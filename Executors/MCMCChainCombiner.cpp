#include<MCMC.hpp>

using namespace std;


int main(int argc, char** argv)
{
	if(argc<3)
	{ 
		cout<<"./MCMCChainCombiner MCMCout.bin MCMC1.bin [MCMC2.bin ... MCMCN.bin]"<<endl;
		return 1;
	}
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();

	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;

	CMCMC MCMC;
	MCMC.CombineBin(argv[1],argv+2, argc-2);
	
	return 0;
}
