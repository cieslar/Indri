#include<MCMC.hpp>
#include<Debug.hpp>

using namespace std;

void PrintHelp()
{
	cout<<"Model builder"<<endl;
	cout<<"-E               - standard pop evo mode"<<endl;
	cout<<"-CAT             - same sim as E-mode but produce artificial catalogue instead of regural output"<<endl;
	cout<<"-T               - tracking mode"<<endl;
	cout<<"-LST             - likelihood stress test mode"<<endl;
	cout<<"-GWPrep          - make a population for the GW processing"<<endl;
	cout<<"Available options:"<<endl;
	cout<<"-n <numofPSR>    - number of simulated pulsars"<<endl;
	cout<<"-i <PopRefFile>  - input file with precomputed geometrical population (issential)"<<endl;
	cout<<"-s <nSeed>       - seed number"<<endl;
	cout<<"-o <Pop.hdf>     - output file (if non specified defualuts to \"Pop.hdf\")"<<endl;
	cout<<"-c <config.xml>  - configuration file"<<endl;
	cout<<"-v               - verbose computations"<<endl;
	cout<<"-h               - prints this help"<<endl;

	cout<<"-k <num>         - number of steps in trackmode or stress mode"<<endl;
	cout<<"-PSRNums <file>  - file with pulsars numbers for tracking"<<endl;
	cout<<"-nNumT <num>     - how many pulsars from PSRNums should be considered"<<endl;
	cout<<"-l <num>         - number of PSR for the Catalogue (output number will not exceed observable PSRs from the model)"<<endl;
	cout<<"-m <Cat.txt>     - the catalogue output"<<endl;

}

int main(int argc, char** argv)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->Sim.bRecycleDynamics=true;
	pCfgPtr->Sim.bDynamics=false;
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bObservations=true;


	if(cmdOptionExists(argv, argv+argc, "-h") || argc < 2) 
	{
		PrintHelp();
		return 1;
	}


	char * input = getCmdOption(argv, argv + argc, "-i");
	if(!input)
	{
		pCfgPtr->Sim.bRecycleDynamics=false;
		pCfgPtr->Sim.bDynamics=true;
	}
	
	//cerr<<pCfgPtr->Sim.bRecycleDynamics<<endl;

	if(cmdOptionExists(argv, argv+argc, "-v"))
	{
		pCfgPtr->Sim.bVerbose = true;
	}
	else
	{
		pCfgPtr->Sim.bVerbose = false;
	}



	char * filename = getCmdOption(argv, argv + argc, "-c");
	if(filename)
	{
		pCfgPtr->ReadXMLMCMC(filename);
	}
	else
	{
		cout<<"Warning: Default setup used"<<endl;
	}




	if(cmdOptionExists(argv, argv+argc, "-s"))
	{
		char * cSeed = getCmdOption(argv, argv + argc, "-s");
		pCfgPtr->Reseed(atoi(cSeed));
	}

	char * output = getCmdOption(argv, argv + argc, "-o");
	if(!output)
	{
		output = new char[256];
		strcpy( output, "Pop.hdf" );
		cout<<"Warning: No output file specified. \"Pop.hdf\" used."<<endl;
	}

	pCfgPtr->PostConfigInit();
	
	if(cmdOptionExists(argv, argv+argc, "-n"))
	{
		char * cNum = getCmdOption(argv, argv + argc, "-n");
		pCfgPtr->Sim.nMaxAvailablePulsarNum=pCfgPtr->Sim.nMaxPulsarNum=atoi(cNum);
	}

    if(	pCfgPtr->Sim.bVerbose )
    {
    	cout<<	pCfgPtr->Sim.nMaxAvailablePulsarNum<< " "<<pCfgPtr->Sim.nMaxPulsarNum<<" "<<pCfgPtr->Sim.nNumberOfIteration<<endl;
    }

	if(cmdOptionExists(argv, argv+argc, "-GWPrep"))
	{
		CModel Model;
		pCfgPtr->Sim.bRecycleDynamics=false;
		pCfgPtr->Sim.bDynamics=true;
		pCfgPtr->Sim.bPhysics=true;
		pCfgPtr->Sim.bObservations=false;
		pCfgPtr->Sim.nNumberOfIteration=1;
		pCfgPtr->Sim.fMaxPulsarAge=50.;
		pCfgPtr->Sim.fMaxDist=250.;
		Model.FullCycle(output);	
	}

	if(cmdOptionExists(argv, argv+argc, "-LST"))
	{
		char * csteps = getCmdOption(argv, argv + argc, "-k");
		if(!csteps)
		{
			cout<<"Error: Number of steps not specified (-k)"<<endl;
			PrintHelp();
			return 3;
		}
		int nsteps(atoi(csteps));

		CModel Model;
		if(input) 
		{
			pCfgPtr->Sim.bRecycleDynamics=true;
			Model.ReadBinary(input);
		}

		ofstream out(output);	
		for(int i(0);i<nsteps;i++)
		{
			Model.FullCycle();
			out<<pCfgPtr->Model.fLnLikelihood<<endl;
		}
		out.close();
	}


	if(cmdOptionExists(argv, argv+argc, "-E"))
	{
		CModel Model;
		if(input) 
		{
			pCfgPtr->Sim.bRecycleDynamics=true;
			Model.ReadBinary(input);
		}

		Model.FullCycle(output);
        if( pCfgPtr->Sim.bVerbose ) cout<<"LnLikelihood: "<<pCfgPtr->Model.fLnLikelihood<<endl;
    
	}

	if(cmdOptionExists(argv, argv+argc, "-T"))
	{
		CPopulation Pop;
		if(pCfgPtr->Sim.bRecycleDynamics) Pop.ReadBinary(input);
		char * ctracks = getCmdOption(argv, argv + argc, "-k");
		if(!ctracks)
		{
			cout<<"Error: Number of steps per track not specified (-k)"<<endl;
			PrintHelp();
			return 3;
		}
		int ntracks=atoi(ctracks);


		char * cnumnums = getCmdOption(argv, argv + argc, "-nNumT");
		if(!cnumnums)
		{
			cout<<"Error: Number of considered pulsars from the PSRNums file not specified"<<endl;
			PrintHelp();
			return 3;
		}
		int nNumNum=atoi(cnumnums);

		int *pNums = new int[nNumNum];

		char * nums = getCmdOption(argv, argv + argc, "-PSRNums");
		if(!nums)
		{
			cout<<"Error: File for PSRs numbers not specified"<<endl;
			PrintHelp();
			return 3;
		}
		ifstream in(nums);
		for(int i(0);i<nNumNum;i++) in>>pNums[i];
		in.close();

		Pop.WriteEvolutionTracksHDF5(output,ntracks,pNums,nNumNum);
		delete [] pNums;
	}


	if(cmdOptionExists(argv, argv+argc, "-CAT"))
	{
		int nCatNum(1000);
		if(cmdOptionExists(argv, argv+argc, "-l"))
		{
			char * cNum = getCmdOption(argv, argv + argc, "-l");
			nCatNum = atoi(cNum);
		}
		else
		{
			cout<<"Warning: No size of the output catalogue provided. Assumed 1000."<<endl;
		}
		
		char * cCatOut = getCmdOption(argv, argv + argc, "-m");
		if(!cCatOut)
		{
			cCatOut = new char[256];
			strcpy( cCatOut, "Catalogue.txt" );
			cout<<"Warning: No catalogue output file specified. \"Catalogue.txt\" used."<<endl;
		}

		CModel Model;
		if(input) 
		{
			pCfgPtr->Sim.bRecycleDynamics=true;
			Model.ReadBinary(input);
		}

		Model.FullCycle();	
		Model.WriteAsCatalogue(cCatOut, nCatNum);
	}

	if(cmdOptionExists(argv, argv+argc, "-LTest"))
    {
        CModel Model;
        CPulsarTester PSRTest;
        PSRTest.CatModelLTest();
    }
	return 0;
}
