#include<MCMC.hpp>
#include<Debug.hpp>
#include<sstream>

using namespace std;

void PrintHelp()
{
	cout<<"MCMCTester"<<endl;
	cout<<"Available options:"<<endl;
	cout<<"-h               - prints this help and provides Example.MCMCConfig.xml"<<endl;
	cout<<"-XMLIO           - XML I/O test"<<endl;
	cout<<"-XMLFishGen      - The xml-o-fish set generator"<<endl;
	cout<<"-HDF             - The hdf5 class test"<<endl;
	cout<<"Test options:"<<endl;
	cout<<"XML I/O:"<<endl;
	cout<<"-c <config.xml>  - input configuration file"<<endl;
	cout<<"-o <config.xml>  - outpu configuration file"<<endl;
	cout<<"-min <num>"<<endl;
	cout<<"-max <num>"<<endl;
	cout<<"-step <num>"<<endl;

}

#include <iostream>
//#include <hdf5.h>
#include <H5Cpp.h>

// Constants
const char saveFilePath[] = "test.h5";
const hsize_t ndims = 1;
const hsize_t ncols = 1024;

void HDFTest()
{
    // Create a hdf5 file
    hid_t file = H5Fcreate(saveFilePath, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    std::cout << "- File created" << std::endl;

    // Create a 2D dataspace
    // The size of the first dimension is unlimited (initially 0)
    // The size of the second dimension is fixed
    hsize_t dims[ndims] = {ncols};
    //hsize_t max_dims[ndims] = {H5S_UNLIMITED};
    hsize_t max_dims[ndims] = {2*ncols};
    hid_t file_space = H5Screate_simple(ndims, dims, max_dims);
    std::cout << "- Dataspace created" << std::endl;

    // Create a dataset creation property list
    // The layout of the dataset have to be chunked when using unlimited dimensions
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    // The choice of the chunk size affects performances
    // This is a toy example so we will choose one line
    hsize_t chunk_dims[ndims] = {ncols};
    H5Pset_chunk(plist, ndims, chunk_dims);
    std::cout << "- Property list created" << std::endl;

    // Create the dataset 'dset1'
    hid_t dset = H5Dcreate(file, "dset1", H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    std::cout << "- Dataset 'dset1' created" << std::endl;

    // Close resources
    H5Pclose(plist);
    H5Sclose(file_space);


    // Create 2D buffer (contigous in memory, row major order)
    // We will allocate enough memory to store 3 lines, so we can reuse the buffer
    float *buffer = new float[ncols];

    for (hsize_t i = 0; i < ncols; ++i) buffer[i]=i;

    // Initial values in buffer to be written in the dataset
    // Create a memory dataspace to indicate the size of our buffer
    // Remember the first buffer is only two lines long
    dims[0] = ncols;
    hid_t mem_space = H5Screate_simple(ndims, dims, NULL);
    std::cout << "- Memory dataspace created" << std::endl;

    // Extend dataset
    // We set the initial size of the dataset to 0x3, we thus need to extend it first
    // Note that we extend the dataset itself, not its dataspace
    // Remember the first buffer is only two lines long
    dims[0] = 2*ncols;
    H5Dset_extent(dset, dims);
    std::cout << "- Dataset extended" << std::endl;



    // Select hyperslab on file dataset
    file_space = H5Dget_space(dset);
    hsize_t start[1] = {0};
    hsize_t count[1] = {ncols};
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    std::cout << "- First hyperslab selected" << std::endl;

    // Write buffer to dataset
    // mem_space and file_space should now have the same number of elements selected
    // note that buffer and &b[0][0] are equivalent
    H5Dwrite(dset, H5T_NATIVE_FLOAT, mem_space, file_space, H5P_DEFAULT, buffer);
    std::cout << "- First buffer written" << std::endl;
    
    // We can now close the file dataspace
    // We could close the memory dataspace and create a new one,
    // but we will simply update its size
    H5Sclose(file_space);

    // ## Second buffer

    // New values in buffer to be appended to the dataset
    for (hsize_t i = 0; i < ncols; ++i) buffer[i]=i*3;
    

    // Resize the memory dataspace to indicate the new size of our buffer
    // The second buffer is three lines long
    //dims[0] = ncols;
    //H5Sset_extent_simple(mem_space, ndims, dims, NULL);
    //std::cout << "- Memory dataspace resized" << std::endl;

    // Extend dataset
    // Note that in this simple example, we know that 2 + 3 = 5
    // In general, you could read the current extent from the file dataspace
    // and add the desired number of lines to it
    //dims[0] = 2*ncols;
    //H5Dset_extent(dset, dims);
    //std::cout << "- Dataset extended" << std::endl;

    // Select hyperslab on file dataset
    // Again in this simple example, we know that 0 + 2 = 2
    // In general, you could read the current extent from the file dataspace
    // The second buffer is three lines long
    file_space = H5Dget_space(dset);
    start[0] = ncols;
    count[0] = ncols;
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    std::cout << "- Second hyperslab selected" << std::endl;

    // Append buffer to dataset
    H5Dwrite(dset, H5T_NATIVE_FLOAT, mem_space, file_space, H5P_DEFAULT, buffer);
    std::cout << "- Second buffer written" << std::endl;
    // Close resources
    delete[] buffer;
    H5Sclose(file_space);
    H5Sclose(mem_space);
    H5Dclose(dset);
    H5Fclose(file);
    std::cout << "- Resources released" << std::endl;
}









void XMLFishGen(int argc, char** argv)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();


	char * InFilename = getCmdOption(argv, argv + argc, "-c");
	if(InFilename)
	{
		pCfgPtr->ReadXMLMCMC(InFilename);
	}
	else
	{
		cout<<"Warning: Default MCMC setup used"<<endl;
	}

	char *cMaxScale = getCmdOption(argv, argv + argc, "-max");
	int nMaxScale;
	if(cMaxScale)
	{
		nMaxScale = atoi(cMaxScale);
	}
	else
	{
		ERROR("Needs the -max argument");
	}


	char *cMinScale = getCmdOption(argv, argv + argc, "-min");
	int nMinScale;
	if(cMinScale)
	{
		nMinScale = atoi(cMinScale);
	}
	else
	{
		ERROR("Needs the -min argument");
	}


	char *cStep = getCmdOption(argv, argv + argc, "-step");
	int nStep;
	if(cStep)
	{
		nStep = atoi(cStep);
	}
	else
	{
		ERROR("Needs the -step argument");
	}




	RealType fBDecayDeltaRangeLog = log10(pCfgPtr->ModelConstraints.Max.fBDecayDelta) - log10(pCfgPtr->ModelConstraints.Min.fBDecayDelta);
	RealType fBDecayKappaRange = pCfgPtr->ModelConstraints.Max.fBDecayKappa - pCfgPtr->ModelConstraints.Min.fBDecayKappa;
	RealType fLPowerLawGammaRangeLog = log10(pCfgPtr->ModelConstraints.Max.fLPowerLawGamma) - log10(pCfgPtr->ModelConstraints.Min.fLPowerLawGamma);
	RealType fLPowerLawAlphaRange = pCfgPtr->ModelConstraints.Max.fLPowerLawAlpha - pCfgPtr->ModelConstraints.Min.fLPowerLawAlpha;
	RealType fLPowerLawBetaRange = pCfgPtr->ModelConstraints.Max.fLPowerLawBeta - pCfgPtr->ModelConstraints.Min.fLPowerLawBeta;
	RealType fPInitMeanRange = pCfgPtr->ModelConstraints.Max.fPInitMean - pCfgPtr->ModelConstraints.Min.fPInitMean;
	RealType fPInitSigmaRange = pCfgPtr->ModelConstraints.Max.fPInitSigma - pCfgPtr->ModelConstraints.Min.fPInitSigma;
	RealType fBInitMeanRange = pCfgPtr->ModelConstraints.Max.fBInitMean - pCfgPtr->ModelConstraints.Min.fBInitMean;
	RealType fBInitSigmaRange = pCfgPtr->ModelConstraints.Max.fBInitSigma - pCfgPtr->ModelConstraints.Min.fBInitSigma;
	RealType fBFinalMeanRange = pCfgPtr->ModelConstraints.Max.fBFinalMean - pCfgPtr->ModelConstraints.Min.fBFinalMean;
	RealType fBFinalSigmaRange = pCfgPtr->ModelConstraints.Max.fBFinalSigma - pCfgPtr->ModelConstraints.Min.fBFinalSigma;
	RealType fDeathPsiRange = pCfgPtr->ModelConstraints.Max.fDeathPsi - pCfgPtr->ModelConstraints.Min.fDeathPsi;



stringstream cIter;

	for(int i(nMinScale);i<=nMaxScale;i+=nStep) if(i>0)
	{
		pCfgPtr->MCMCStepSize.fBDecayDelta = pow(10.,fBDecayDeltaRangeLog/static_cast<float>(i));
		pCfgPtr->MCMCStepSize.fBDecayKappa = fBDecayKappaRange/static_cast<float>(i);
		pCfgPtr->MCMCStepSize.fLPowerLawGamma = pow(10.,fLPowerLawGammaRangeLog/static_cast<float>(i));
		pCfgPtr->MCMCStepSize.fLPowerLawAlpha = fLPowerLawAlphaRange/static_cast<float>(i);
		pCfgPtr->MCMCStepSize.fLPowerLawBeta = fLPowerLawBetaRange/static_cast<float>(i);
		pCfgPtr->MCMCStepSize.fBInitMean = fBInitMeanRange/static_cast<float>(i);
		pCfgPtr->MCMCStepSize.fBInitSigma = fBInitSigmaRange/static_cast<float>(i);
		pCfgPtr->MCMCStepSize.fPInitMean = fPInitMeanRange/static_cast<float>(i);
		pCfgPtr->MCMCStepSize.fPInitSigma = fPInitSigmaRange/static_cast<float>(i);
		pCfgPtr->MCMCStepSize.fDeathPsi = fDeathPsiRange/static_cast<float>(i);

		cIter.str("");
		cIter.clear();
		cIter<<"MCMC."<<i<<"nthFishpart.xml";
		pCfgPtr->WriteXMLMCMC(cIter.str().c_str());
	}







}



void XMLInOutTest(int argc, char** argv)
{

	CConfiguration *pCfgPtr = CConfiguration::GetInstance();


	char * InFilename = getCmdOption(argv, argv + argc, "-c");
	if(InFilename)
	{
		pCfgPtr->ReadXMLMCMC(InFilename);
	}
	else
	{
		cout<<"Warning: Default MCMC setup used"<<endl;
	}

	char * OutFilename = getCmdOption(argv, argv + argc, "-o");
	if(OutFilename)
	{
		pCfgPtr->WriteXMLMCMC(OutFilename);
	}
	else
	{
		pCfgPtr->WriteXMLMCMC("Output.MCMCConfig.xml");
	}


}
int main(int argc, char** argv)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	if(cmdOptionExists(argv, argv+argc, "-h") || argc==1) 
	{
		PrintHelp();
		pCfgPtr->WriteXMLMCMC("Example.MCMCConfig.xml");
		return 1;
	}

	if(cmdOptionExists(argv, argv+argc, "-XMLIO")) 
	{
		XMLInOutTest(argc,argv);
	}


	if(cmdOptionExists(argv, argv+argc, "-XMLFishGen")) 
	{
		XMLFishGen(argc,argv);
	}

	if(cmdOptionExists(argv, argv+argc, "-HDF")) 
	{
		HDFTest();
	}

	return 0;
}
