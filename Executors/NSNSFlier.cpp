#include<Galaxy.hpp>
#include<PulsarDynamics.hpp>
#include<iomanip>
#include<fstream>
#include<string>
#include<ConfigContainer.hpp>
#include<Statistics.hpp>
#include<MCMC.hpp>
#include<NSNS.hpp>

#include<Population.hpp>
using namespace std;



void PrintHelp()
{
	cout<<"NSNSFlier flies you through the Galaxy like a NS-NS binary."<<endl;
	cout<<"./NSNSFlier -i InputFile -o OutputFile [-s nSeed]"<<endl<<endl;
	cout<<"-i InputFile - just a filename. Essential."<<endl;
	cout<<"-o OutpuFile - just a filename. Essential."<<endl;
	cout<<"-s nSeed - Seed for the random generator. Optional."<<endl;
	cout<<"Structer of Input File:"<<endl;
	cout<<"nIdNumber fTBoring[Myr] fT1[Myr] fT2[Myr] fV1KickX[km/s] fV1KickY[km/s] fV1KickZ[km/s] fV2KickX[km/s] fV2KickY[km/s] fV2KickZ[km/s]"<<endl<<endl;
	cout<<"fTBoring - the evolution time from the beggining to the first kick. Extremly boring..."<<endl;
	cout<<"fT1 - the evolution time between the first and the second kick."<<endl;
	cout<<"fT2 - the evolution time between the second kick and the observation."<<endl;
	cout<<"Programe assumes that each lines ends with standard unix EoL."<<endl;
}



int main(int argc, char** argv)
{
	if(argc==1)
	{
		PrintHelp();
		return 1;
	}
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit(); 
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bPhysics=false;
	pCfgPtr->Sim.bObservations=false;
	pCfgPtr->Sim.bVerbose=false;    

	char *InputFilename = getCmdOption(argv, argv + argc, "-i");
	char *OutputFilename = getCmdOption(argv, argv + argc, "-o");

	if(cmdOptionExists(argv, argv+argc, "-s"))
	{
		char *cSeed = getCmdOption(argv, argv + argc, "-s");
		pCfgPtr->Reseed(atoi(cSeed));
	}



	CNSNSBinaryPopulation NSNS;
	NSNS.Read(InputFilename);

	NSNS.Init();
	NSNS.Evolve();

	NSNS.Project();

	NSNS.Write(OutputFilename);
	NSNS.WriteProjectionMtx();

	return 0;
}
