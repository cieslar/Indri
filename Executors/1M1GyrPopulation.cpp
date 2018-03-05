#include<iostream>
#include<Population.hpp>
#include<Debug.hpp>


using namespace std;

int main(int argc, char ** argv)
{

#if 0
	if(argc1)
	{
		cout<<"Available options:"<<endl;
		cout<<"-s <number>     Custom seed for the random engine (positive intiger)."<<endl;
		cout<<"-f <filename>   Custom filename."<<endl;
	}
#endif
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();


	pCfgPtr->Sim.nMaxPulsarNum=1000000;//1M
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=true;
	pCfgPtr->Sim.fMaxPulsarAge=20.;//1Gyr


	char * filename = getCmdOption(argv, argv + argc, "-f");
	if(!filename)
	{
		filename = new char[256];
		strcpy( filename, "1M1GyrPopulation.dat" );
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
	Pop.Init();

        Pop.Evolve();


	Pop.WriteBinary(filename);


	return 0;
}


