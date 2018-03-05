#include<MCMC.hpp>

using namespace std;


int main(int argc, char** argv)
{
	if(argc!=3)
	{
		cout<<"./MCMC MCMC.bin MCMC.hdf"<<endl;
		return 1;
	}
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();


	CMCMC MCMC;
	cout<<"Reading"<<endl;
	MCMC.ReadBinary(argv[1]);

	cout<<"Writing"<<endl;
	MCMC.WriteBinaryHDF5(argv[2]);


	return 0;
}
