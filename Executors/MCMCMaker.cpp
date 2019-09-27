#include<MCMC.hpp>
#include<Debug.hpp>

using namespace std;

void PrintHelp()
{
	cout<<"Monte chain Markov Carlo for Pulsar Population"<<endl;
	cout<<"Available options:"<<endl;
	cout<<"-i <PopRefFile>  - input file with precomputed geometrical population (issential)"<<endl;
	cout<<"-s <nSeed>       - seed number"<<endl;
	cout<<"-o <MCMCFile>    - output file (if non specified defualuts to \"MCMCOut.bin\")"<<endl;
	cout<<"-c <config.xml>  - configuration file"<<endl;
	cout<<"-v               - verbose MCMC computations"<<endl;
	cout<<"-vv              - verbose MCMC and Model computations"<<endl;
	cout<<"-h               - prints this help and provides Example.MCMCConfig.xml"<<endl;
    cout<<"-n               - number of pulsars to consider"<<endl;
    cout<<"-m <MCMCFile>    - the .bin MCMC file, in case of continuing from a given point in the chain"<<endl;
    cout<<"-k <LinkNum>     - the point in chain (array-wise) to continoue from"<<endl;
    cout<<"-l <ChainNum>    - the MCMC chain for continuation"<<endl;
    cout<<"-r               - switches on the random param walker"<<endl;
}

int main(int argc, char** argv)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->Sim.nNumOfChains=1;
	pCfgPtr->Sim.nNumPerChain=10;

	if(cmdOptionExists(argv, argv+argc, "-h")) 
	{
		PrintHelp();
		pCfgPtr->WriteXMLMCMC("Example.MCMCConfig.xml");
		return 1;
	}

	pCfgPtr->Sim.bMCMCVerbose = false;
	pCfgPtr->Sim.bVerbose = false;
	if(cmdOptionExists(argv, argv+argc, "-v"))
	{
		pCfgPtr->Sim.bMCMCVerbose = true;
	}

	if(cmdOptionExists(argv, argv+argc, "-vv"))
	{
		pCfgPtr->Sim.bMCMCVerbose = true;
		pCfgPtr->Sim.bVerbose = true;
	}


    char * filename = getCmdOption(argv, argv + argc, "-c");
	if(filename)
	{
		pCfgPtr->ReadXMLMCMC(filename,pCfgPtr->Sim.bMCMCVerbose);
	}
	else
	{
		cout<<"Warning: Default MCMC setup used"<<endl;
	}
	pCfgPtr->Sim.bRecycleDynamics=true;
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bObservations=true;








	char * input = getCmdOption(argv, argv + argc, "-i");
	if(!input)
	{
		cout<<"Error: No input file specified"<<endl;
		PrintHelp();
		return 2;
	}








	if(cmdOptionExists(argv, argv+argc, "-s"))
	{
		char * cSeed = getCmdOption(argv, argv + argc, "-s");
		pCfgPtr->Reseed(atoi(cSeed));
	}




	if(cmdOptionExists(argv, argv+argc, "-Seq"))
	{
		pCfgPtr->Sim.bMCMCSequential = true;
	}

	if(cmdOptionExists(argv, argv+argc, "-n"))
	{
		char * cNum = getCmdOption(argv, argv + argc, "-n");
		pCfgPtr->Sim.nMaxAvailablePulsarNum=pCfgPtr->Sim.nMaxPulsarNum=atoi(cNum);
	}


    

	if(cmdOptionExists(argv, argv+argc, "-r"))
	{
		pCfgPtr->Sim.bUseRandomParamSampler=true;
    }





	pCfgPtr->PostConfigInit();
	

	
	CMCMC MCMC;



    if(cmdOptionExists(argv, argv+argc, "-m"))
    {
        pCfgPtr->Sim.bContinueMCMCFromLastPostion=true;
        char *cMCMCChainsInput = getCmdOption(argv, argv + argc, "-m");

        if(cmdOptionExists(argv, argv+argc, "-k"))
        {
            char *cLastMCMCPostion = getCmdOption(argv, argv + argc, "-k");
            pCfgPtr->Sim.nLastMCMCPostion = atoi(cLastMCMCPostion);
        }
        else
        {
            pCfgPtr->Sim.nLastMCMCPostion=0;
            cerr<<"The last MCMC position not specified. Using: "<<pCfgPtr->Sim.nLastMCMCPostion<<endl;
        }

        if(cmdOptionExists(argv, argv+argc, "-l"))
        {
            char *cLastMCMCChain = getCmdOption(argv, argv + argc, "-l");
            pCfgPtr->Sim.nLastMCMCChain = atoi(cLastMCMCChain);
        }
        else
        {
            pCfgPtr->Sim.nLastMCMCChain=0;
            cerr<<"The MCMC chain for continuation not specified not specified. Using: "<<pCfgPtr->Sim.nLastMCMCChain<<endl;
        }

        MCMC.ReadBinary(cMCMCChainsInput);
    }


	MCMC.PopReadBinary(input);
	if(pCfgPtr->Sim.bMCMCSequential)
	{
		////////////////////
		//Work in progress//
		////////////////////
	}
	else
	{
		MCMC.MakeChains();
	}
	char * output = getCmdOption(argv, argv + argc, "-o");
	if(!output)
	{
		output = new char[256];
		strcpy( output, "MCMCOut.bin" );
		cout<<"Warning: No output file specified. \"MCMCOut.bin\" used."<<endl;
	}

	MCMC.WriteBinary(output);

	return 0;
}
