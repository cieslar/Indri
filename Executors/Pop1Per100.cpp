#include<iostream>
#include<Population.hpp>

#include <algorithm>
#include<string>
#include<cstdlib>

using namespace std;

void PrintHelp()
{
	cout<<"========== Pop1Per100 v2 =========="<<endl;
	cout<<"Fly your Neutron Stars through Galactic potential today! Now with Hobbs' kicks!"<<endl;
	cout<<"Population of evenly spaced (1 per 100 years) pulsars that propagate through three-component Galactic potential. Optional P-Pdot-B evolution and Parkes Multibeam survey."<<endl;
	cout<<"Available options:"<<endl;
	cout<<"-n <number>          Number of simulated pulsars. Essential option."<<endl;
	cout<<"-v                   Verbose console output."<<endl;
	cout<<"-s <number>          Custom seed for the random engine (positive intiger)."<<endl;
	cout<<"-f <filename>        Custom output filename."<<endl;
	cout<<"-SigAge <sig>        Smear the age using Gauss with sigma sig [years]."<<endl;
	cout<<"-PhysEvo             Switch on the P-Pdot-B evolution."<<endl;
	cout<<"-xml <config.xml>    Xml config file. Relevant option will be overwritten by option flags."<<endl;
	cout<<"-Obs                 Switch on the ParkesMb survey."<<endl;
	cout<<"-h                   Prints this help and provides Example.MCMCConfig.xml"<<endl;

}



int main(int argc, char ** argv)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();

	if(argc<2)
	{
		PrintHelp();
		return 1;
	}

	if(cmdOptionExists(argv, argv+argc, "-h")) 
	{
		PrintHelp();
		pCfgPtr->WriteXMLMCMC("Example.MCMCConfig.xml");
		return 1;
	}


	char * xmlfilename = getCmdOption(argv, argv + argc, "-xml");
	if(xmlfilename)
	{
		pCfgPtr->ReadXMLMCMC(xmlfilename);
	}
	else
	{
		cout<<"No xml specified. Default config used."<<endl;
	}


	if(cmdOptionExists(argv, argv+argc, "-PhysEvo"))
	{
		pCfgPtr->Sim.bPhysics=true;
	}
	else
	{
		pCfgPtr->Sim.bPhysics=false;
	}


	if(cmdOptionExists(argv, argv+argc, "-Obs"))
	{
		pCfgPtr->Sim.bObservations=true;
	}
	else
	{
		pCfgPtr->Sim.bObservations=false;
	}






	if(cmdOptionExists(argv, argv+argc, "-v"))
	{
		pCfgPtr->Sim.bVerbose = true;
	}
	else
	{
		pCfgPtr->Sim.bVerbose = false;
	}

	char * filename = getCmdOption(argv, argv + argc, "-f");
	if(!filename)
	{
		filename = new char[16];
		strcpy( filename, "Pop1Per100.txt" );
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
		pCfgPtr->Sim.nMaxAvailablePulsarNum=atoi(cNum);
	}
	else
	{
		cerr<<"Error! Number of simulated pulsars not specified"<<endl;
		PrintHelp();
		return 2;
	}



 	if(pCfgPtr->Sim.bVerbose) cout<<"Population init."<<endl;
	CPopulation Pop;

	char * cSigma = getCmdOption(argv, argv+argc, "-SigAge");
	if(cSigma)
	{
		Pop.Init1Per100(true,atof(cSigma));
	}
	else Pop.Init1Per100();


 	if(pCfgPtr->Sim.bVerbose) cout<<"Evolving"<<endl;

        Pop.Evolve();

 	if(pCfgPtr->Sim.bVerbose) cout<<"Writing results to "<<filename<<endl;
	
	Pop.Write(filename);

 	
	return 0;
}


