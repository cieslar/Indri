#include<iostream>
#include<Population.hpp>

using namespace std;


//cuts the population


void PrintHelp()
{
	cout<<"Various postprocessing for population output"<<endl;
	cout<<"Major modes:"<<endl;
	cout<<"-C    cut a subset of given max size and survey geometry"<<endl;
	cout<<"-T    translate binary population's output to hdf5 file"<<endl;
	cout<<"-M    modify population (currently hardcoded)"<<endl;
	cout<<"-h    print this help"<<endl;
	cout<<endl;
	cout<<"Options:"<<endl;
	cout<<"-NOut <number> maximum population size"<<endl;
	cout<<"-i <filename>  input file"<<endl;
	cout<<"-o <filename>  output file"<<endl;
	cout<<"TODO:"<<endl;
	cout<<"modification as a part of the xml"<<endl;
}

int main(int argc, char** argv)
{
	if((argc==1) || cmdOptionExists(argv, argv+argc, "-h"))
	{
		PrintHelp();
		return 0;
	}
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	char * ile = getCmdOption(argv, argv + argc, "-NOut");
	if(!ile)	
	{
		cout<<"Warning: Nax population size not specified"<<endl;
		cout<<"Using nMaxAvailablePulsarNum: "<<pCfgPtr->Sim.nMaxAvailablePulsarNum<<endl;
		pCfgPtr->Sim.nMaxAvailablePulsarNum=pCfgPtr->Sim.nMaxPulsarNum;
	}
	else
	{
		pCfgPtr->Sim.nMaxAvailablePulsarNum=pCfgPtr->Sim.nMaxPulsarNum=atoi(ile);
	}

	char * input = getCmdOption(argv, argv + argc, "-i");
	if(!input)
	{
		cout<<"Error: No input file specified"<<endl;
		PrintHelp();
		return 2;
	}


	char * output = getCmdOption(argv, argv + argc, "-o");
	if(!output)
	{
		cout<<"Error: No output file specified"<<endl;
		PrintHelp();
		return 3;
	}
	
	pCfgPtr->PostConfigInit();

	CPopulation Pop;

	if(cmdOptionExists(argv, argv+argc, "-M"))
	{
		Pop.ReadBinaryRelevantGeometry(input);
		for(int i(pCfgPtr->Sim.nMaxPulsarNum);i--;) Pop.pPulsar[i].fAge/=2.5;
		Pop.WriteBinary(output);
	}


	if(cmdOptionExists(argv, argv+argc, "-C")) 
	{
		Pop.ReadBinaryRelevantGeometry(input);
		Pop.WriteBinary(output);
	}

	if(cmdOptionExists(argv, argv+argc, "-T")) 
	{
		Pop.ReadBinary(input);
		Pop.WriteHDF5(output,0);
	}

	return 0;
}
