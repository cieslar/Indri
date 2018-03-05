#include<MCMC.hpp>
#include<Debug.hpp>

using namespace std;

void PrintHelp()
{
	cout<<"DM Changer"<<endl;
	cout<<"Modes:"<<endl;
	cout<<"-prep                  - read the popref and preapare a DMfile for the YMW16 computation (exteran script)"<<endl;
	cout<<"-change                - substitues the DM in the popref according to the YMW16 DM-file."<<endl;
	cout<<"Available options:"<<endl;
	cout<<"-i <PopRefFile>        - input file with precomputed geometrical population (issential)"<<endl;
	cout<<"-o <ModPopRefFile>     - output file"<<endl;
	cout<<"-idm <YMW16DMFile>     - input DM file"<<endl;
	cout<<"-odm <DMFile>          - output DM file"<<endl;
	cout<<"-h                     - print this help"<<endl;
	cout<<"Usage:"<<endl;
	cout<<"-prep -i <PopRefFile> -odm <DMFile>"<<endl;
	cout<<"-change -i <PopRefFile> -idm <YMW16DMFile> -o <ModPopRefFile>"<<endl;
	

}

int main(int argc, char** argv)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();

	if(cmdOptionExists(argv, argv+argc, "-h") || argc < 2) 
	{
		PrintHelp();
		return 1;
	}



	char * PopRefFile = getCmdOption(argv, argv + argc, "-i");
	if(!PopRefFile) ERROR("");


//	pCfgPtr->PostConfigInit();

	CPopulation Pop;
	Pop.ReadBinary(PopRefFile);		

	if(cmdOptionExists(argv, argv+argc, "-prep"))
	{
		char * ODMFile = getCmdOption(argv, argv + argc, "-odm");
		if(!ODMFile) ERROR("");
		Pop.WriteTxtDM(ODMFile);
	}


	if(cmdOptionExists(argv, argv+argc, "-change"))
	{
		char * IDMFile = getCmdOption(argv, argv + argc, "-idm");
		if(!IDMFile) ERROR("");

		char * ModPopRefFile = getCmdOption(argv, argv + argc, "-o");
		if(!ModPopRefFile) ERROR("");

		Pop.ReadTxtDM(IDMFile);
		Pop.WriteBinary(ModPopRefFile);
	}



	return 0;
}
