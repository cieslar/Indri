#include<MCMC.hpp>
#include<Debug.hpp>

using namespace std;

void PrintHelp()
{
	cout<<"A testing tool for the xml configuration."<<endl;
	cout<<"Available options:"<<endl;
	cout<<"-oc [<MCMCFile>] - writes the xml file (with optional filename, default set to MCMCConfig.xml)"<<endl;
	cout<<"-ic <config.xml> - reads the configuration file"<<endl;
    cout<<"-s <scale>       - scales the parameters' jumps with the value scale"<<endl;
    cout<<"-h               - prints this help"<<endl;
}

int main(int argc, char** argv)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();

	if( (argc<2) || cmdOptionExists(argv, argv+argc, "-h") ) 
	{
		PrintHelp();
		return 1;
	}


	if(cmdOptionExists(argv, argv+argc, "-ic"))
    { 
        char * filename = getCmdOption(argv, argv + argc, "-ic");
	    if(filename)
	    {
		    pCfgPtr->ReadXMLMCMC(filename);
	    }
	    else
	    {
		    cout<<"Error: No filename specified."<<endl;
            return 2;
    	}
    }
    if(cmdOptionExists(argv,argv+argc, "-s"))
    {
        char * cScale = getCmdOption(argv, argv + argc, "-s");
        RealType fScale = atof(cScale);
        pCfgPtr->Sim.fJumpScallingFactor=fScale;
        CModelParameterDistributions Params;
        Params.ScaleJumpSigmas();
    }
 	if(cmdOptionExists(argv, argv+argc, "-oc"))
    { 
        char * outfilename = getCmdOption(argv, argv + argc, "-oc");
	    if(outfilename)
	    {
		    pCfgPtr->WriteXMLMCMC(outfilename);
	    }
	    else
	    {
		    cout<<"Warning: No xml output filename specified. Default MCMCConfig.xml used."<<endl;
            pCfgPtr->WriteXMLMCMC("MCMCConfig.xml");
    	}
    }

    return 0;
}
