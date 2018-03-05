#include<MCMC.hpp>
#include<sstream>
#include<iostream>

using namespace std;
int gnModelsPerParam(1);
int gnMaxPulsarNum(10000);
char gcPopRef[256];



void TestMaxwell()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();


	CTestGalaxy Gal;
	Gal.PaczynksiKickTest(gnMaxPulsarNum);
	Gal.HobbsKickTest(gnMaxPulsarNum);


}


void TestDMDiff()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	//pCfgPtr->Sim.bUseDMContainer=true;
	//pCfgPtr->Sim.bUsePrecomputedDMData=true;
	pCfgPtr->Sim.bComputeDMDiff=true;

	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	CModel Model;



	Model.Zero();
	Model.Evolve();

	//Model.WriteBinaryHDF5("DMDiffPop.hdf");


}


void TestDistributions()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();


	CPulsarPhysicTester PSR;

	PSR.InitPeriod(100000,"InitialPeriodDistribution.dat");
	PSR.FinalBField(100000,"MinimalBFieldDistribution.dat");
	PSR.InitBField(100000,"InitialBFieldDistribution.dat");

}

void TestPulsarDeathArea()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();


	CPulsarTester PSR;

	PSR.DeathTester();

}

void TestLErot()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();


	CPulsarTester PSR;



	PSR.LErotTester();
}

void TestATNFCoverage()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=true;
	pCfgPtr->Sim.bUseCatalogue=true;
	pCfgPtr->PostConfigInit();

	CObservations Obs(ParkesMultiBeam);
}

void TestDMContainer()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
pCfgPtr->Sim.bUsePrecomputedDMData=false;
	pCfgPtr->Sim.bVerbose=true;
	CDMContainer *pDM=pCfgPtr->pDM;

	//pDM->WriteBinary();
}


void PrintFilteredATNF()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();

	pCfgPtr->Sim.bVerbose=true;

	CObservations ATNF(ParkesMultiBeam);
	ATNF.ComputeMaxDist();

	

}
void TestMaxDistATNF()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();

	pCfgPtr->Sim.bVerbose=true;

	CObservations ATNF(ParkesMultiBeam);
	ATNF.ComputeMaxDist();

	

}
//A test it the Narayan Ostriker Luminosity makes population visible in Parkes
void TestModelTimeDMContainer()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.bUseDMContainer=true;
	pCfgPtr->Sim.bUsePrecomputedDMData=true;
	
	//temp - aby sprawdzic czy sa konsystentne
	//pCfgPtr->Sim.LuminosityModel=PowerLaw;


	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;
TTimestamp A,B;
RealType Time;
int nModels(0);

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;

	A=GetTimestamp();
	Model.Zero();
	Model.Evolve();
	Model.Densify();
	fLnL=Model.LnLikelihood();
	B=GetTimestamp();

	cout<<B-A<<endl;
}
void TestModelStability()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.bUseDMContainer=true;
	pCfgPtr->Sim.bUsePrecomputedDMData=true;
	
	//temp - aby sprawdzic czy sa konsystentne
	//pCfgPtr->Sim.LuminosityModel=PowerLaw;


	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;
TTimestamp A,B,C;
RealType Time;
int nModels(0);

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;

	A=GetTimestamp();
	Model.Zero();
	Model.Evolve();
	Model.Densify();
	fLnL=Model.LnLikelihood();
	B=GetTimestamp();
	//Model.WriteBinaryHDF5("dupa.hdf",0);
	C=GetTimestamp();

	Model.PrintPulsar(0);	
	cout<<B-A<<" "<<C-B<<endl;
}
void TestModelTimeNE2001()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.bUseDMContainer=false;
	pCfgPtr->Sim.bUsePrecomputedDMData=false;
	
	//temp - aby sprawdzic czy sa konsystentne
	//pCfgPtr->Sim.LuminosityModel=PowerLaw;


	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;
TTimestamp A,B;
RealType Time;
int nModels(0);

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;

	A=GetTimestamp();
	Model.Zero();
	Model.Evolve();
	Model.Densify();
	fLnL=Model.LnLikelihood();
	B=GetTimestamp();

	cout<<B-A<<endl;
}


void TestStackPrint()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();

//	pCfgPtr->ReadXML("GalacticPopulation.xml");
	cout<<pCfgPtr->pRan->Gauss(-0.01,0.)<<endl;
}


void VisibilityTestParkes()
{
	CSurvey Parkes(ParkesMultiBeam);	
	Galactic Gal;
	for(Gal.fL=0.;Gal.fL<=360.;Gal.fL+=1.) for(Gal.fB=-90.;Gal.fB<=90.;Gal.fB+=1.)
	{
		if(Parkes.IsCovered( Gal.fL, Gal.fB)) cout<<Gal.fL<<" "<<Gal.fB<<" "<<Parkes.IsCovered( Gal.fL, Gal.fB)<<endl;
	}

}



void TransformTest(double x, double y, double z)
{
	Galactic Gal;
	Cartesian LocalPos;
	LocalPos.fX=x;
	LocalPos.fY=y;
	LocalPos.fZ=z;

	RealType fDistPlane = sqrt(LocalPos.fX*LocalPos.fX + LocalPos.fY*LocalPos.fY);

	Gal.fDist=Length(LocalPos);
	Gal.fB = atan2(LocalPos.fZ, fDistPlane);

	if(Gal.fB > M_PI/2. || Gal.fB < -M_PI/2.)
	{
		ERROR("Coordinates' transformation is wrong.");
	}

	Gal.fL = atan2(LocalPos.fY,LocalPos.fX);
	if(Gal.fL<0.) Gal.fL+=2.*M_PI;

	cout<<Gal.fL*180./M_PI<<" "<<Gal.fB*180./M_PI<<" "<<Gal.fDist<<endl;
}


void ConfigTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.bMCMCVerbose=true;

	pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	pCfgPtr->MCMCParameterUse.bLPowerLawGamma = false;
	pCfgPtr->MCMCParameterUse.bLPowerLawAlpha = false; 
	pCfgPtr->MCMCParameterUse.bLPowerLawBeta = false;
	pCfgPtr->MCMCParameterUse.bBDecayDelta = true;
	pCfgPtr->MCMCParameterUse.bBDecayKappa = false;

	pCfgPtr->MCMCParameterConstraint.bLPowerLawGamma = false;
	pCfgPtr->MCMCParameterConstraint.bLPowerLawAlpha = false; 
	pCfgPtr->MCMCParameterConstraint.bLPowerLawBeta = false;
	pCfgPtr->MCMCParameterConstraint.bBDecayDelta = true;
	pCfgPtr->MCMCParameterConstraint.bBDecayKappa = false;
	
	pCfgPtr->Sim.nNumOfChains=100;
	pCfgPtr->Sim.nNumPerChain=100;

	pCfgPtr->WriteXML("test.xml");


}

void NarayanLDeathPsiTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;

//	pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	pCfgPtr->MCMCParameterUse.bDeathPsi=true;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fPsi( pCfgPtr->ModelConstraints.Max.fDeathPsi );
	int nIle(gnModelsPerParam);
	while(fPsi >= pCfgPtr->ModelConstraints.Min.fDeathPsi)
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fDeathPsi=fPsi;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fPsi<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fPsi-=.01;
	}

}




void NarayanLBInitMeanTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;

//	pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;
	pCfgPtr->MCMCParameterUse.bBInitMean=true;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fMean( pCfgPtr->ModelConstraints.Max.fBInitMean );
	int nIle(gnModelsPerParam);
	while(fMean >= pCfgPtr->ModelConstraints.Min.fBInitMean)
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fBInitMean=fMean;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fMean<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fMean-=.01;
	}

}




void SeqentionaParameterTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.bRecycleDynamics=true;

//	pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	pCfgPtr->MCMCParameterUse.bBInitSigma=true;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fD,fAvgD,fMinD,fMaxD;
	int nIle(gnModelsPerParam);
	
	RealType fSigma( pCfgPtr->ModelConstraints.Max.fBInitSigma );
	while(fSigma >= pCfgPtr->ModelConstraints.Min.fBInitSigma)
	{
		fAvgLnL=0.;
		fAvgD=0.;


		pCfgPtr->Model.fBInitSigma=fSigma;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		Model.ComputeMinimizationValue();

		fLnL=pCfgPtr->Model.fLnLikelihood;
		fD  =pCfgPtr->Model.fDValue;

		fAvgLnL+=fLnL;
		fAvgD  +=fD;

		fMinLnL=fMaxLnL=fLnL;
		fMinD=fMaxD=fD;

		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			Model.ComputeMinimizationValue();

			fLnL=pCfgPtr->Model.fLnLikelihood;
			fD  =pCfgPtr->Model.fDValue;
			
			fAvgLnL+=fLnL;
			fAvgD  +=fD;

			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMinD>fD)     fMinD  =fD;

			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
			if(fMaxD  <fD)   fMaxD  =fD;
		}
		fAvgLnL/=nIle;
		fAvgD  /=nIle;

		cout<<fSigma<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL;//<<endl;
		cout<<" "<<fAvgD<<" "<<fMinD<<" "<<fMaxD<<endl;
		fSigma-=.01;
	}


}



void NarayanLBInitSigmaTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.bRecycleDynamics=true;

//	pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	pCfgPtr->MCMCParameterUse.bBInitSigma=true;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fSigma( pCfgPtr->ModelConstraints.Max.fBInitSigma );
	int nIle(gnModelsPerParam);
	while(fSigma >= pCfgPtr->ModelConstraints.Min.fBInitSigma)
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fBInitSigma=fSigma;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fSigma<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fSigma-=.01;
	}

}

void NarayanLPInitMeanTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;

//	pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	pCfgPtr->MCMCParameterUse.bPInitMean=true;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fMean( pCfgPtr->ModelConstraints.Max.fPInitMean );
	int nIle(gnModelsPerParam);
	while(fMean >= pCfgPtr->ModelConstraints.Min.fPInitMean)
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fPInitMean=fMean;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fMean<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fMean-=.01;
	}

}

void NarayanLPInitSigmaTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;

//	pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	pCfgPtr->MCMCParameterUse.bPInitSigma=true;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fSigma( pCfgPtr->ModelConstraints.Max.fPInitSigma );
	int nIle(gnModelsPerParam);
	while(fSigma >= pCfgPtr->ModelConstraints.Min.fPInitSigma)
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fPInitSigma=fSigma;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fSigma<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fSigma-=.01;
	}

}
void NarayanLBDeltaMCMCTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=true;
	pCfgPtr->Sim.bMCMCVerbose=true;

	//pCfgPtr->Sim.fMaxPulsarAge=25.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	pCfgPtr->MCMCParameterUse.bLPowerLawGamma = false;
	pCfgPtr->MCMCParameterUse.bLPowerLawAlpha = false; 
	pCfgPtr->MCMCParameterUse.bLPowerLawBeta = false;
	pCfgPtr->MCMCParameterUse.bBDecayDelta = true;
	pCfgPtr->MCMCParameterUse.bBDecayKappa = false;

	pCfgPtr->MCMCParameterConstraint.bLPowerLawGamma = false;
	pCfgPtr->MCMCParameterConstraint.bLPowerLawAlpha = false; 
	pCfgPtr->MCMCParameterConstraint.bLPowerLawBeta = false;
	pCfgPtr->MCMCParameterConstraint.bBDecayDelta = true;
	pCfgPtr->MCMCParameterConstraint.bBDecayKappa = false;
	
	pCfgPtr->Sim.nNumOfChains=100;
	pCfgPtr->Sim.nNumPerChain=100;



	CMCMCPostProcessing MCMCTest;
	MCMCTest.MakeChains();
	MCMCTest.WriteBinary("MCMCTest1p.bin");
	MCMCTest.MakeHistograms();




}


void PowerLawLRotTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.LuminosityModel=PowerLaw;
	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fAlpha=pCfgPtr->ModelConstraints.Max.fLPowerLawAlpha;
	RealType fBeta = -(1./3.)*fAlpha;

	int nIle(gnModelsPerParam);
	while(fAlpha >= pCfgPtr->ModelConstraints.Min.fLPowerLawAlpha )
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fLPowerLawAlpha=fAlpha;
		pCfgPtr->Model.fLPowerLawAlpha=fBeta;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fAlpha<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fAlpha-=.01;
		fBeta = -(1./3.)*fAlpha;
	}

}



void PowerLawLGammaTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.LuminosityModel=PowerLaw;
	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fGamma=pCfgPtr->ModelConstraints.Max.fLPowerLawGamma;

	int nIle(gnModelsPerParam);
	while(fGamma >= pCfgPtr->ModelConstraints.Min.fLPowerLawGamma )
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fLPowerLawGamma=fGamma;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fGamma<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fGamma/=1.05;
	}

}
void PowerLawLAlphaTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.LuminosityModel=PowerLaw;
	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fAlpha=pCfgPtr->ModelConstraints.Max.fLPowerLawAlpha;

	int nIle(gnModelsPerParam);
	while(fAlpha >= pCfgPtr->ModelConstraints.Min.fLPowerLawAlpha )
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fLPowerLawAlpha=fAlpha;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fAlpha<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fAlpha-=.01;
	}

}
void PowerLawLBetaTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.LuminosityModel=PowerLaw;
	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fBeta=pCfgPtr->ModelConstraints.Max.fLPowerLawBeta;

	int nIle(gnModelsPerParam);
	while(fBeta >= pCfgPtr->ModelConstraints.Min.fLPowerLawBeta )
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fLPowerLawBeta=fBeta;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fBeta<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fBeta-=.01;
	}

}

void NarayanLBKappaSeqTest()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;

	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;

	CModel Model;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
//pCfgPtr->ModelConstraints.Max.fBDecayKappa=1.6;
//pCfgPtr->ModelConstraints.Min.fBDecayKappa=-1.;
	RealType fKappa=pCfgPtr->ModelConstraints.Max.fBDecayKappa;

	int nIle(gnModelsPerParam);
//	for(int i(0);i<10;i++)
	while(fKappa >= pCfgPtr->ModelConstraints.Min.fBDecayKappa )
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fBDecayKappa=fKappa;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fKappa<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fKappa-=0.01;
	}


}





//A test it the Narayan Ostriker Luminosity makes population visible in Parkes
void NarayanLBDeltaSeqTest(const char* PopRef=NULL)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.bMCMCVerbose=false;
	pCfgPtr->Sim.bRecycleDynamics=true;



	//temp - aby sprawdzic czy sa konsystentne
	//pCfgPtr->Sim.LuminosityModel=PowerLaw;


	//pCfgPtr->Sim.fMaxPulsarAge=5.;
	//pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum;
	
	CModel Model;

	if(PopRef) 
	{
		pCfgPtr->Sim.bDynamics=false;
		pCfgPtr->Sim.bRecycleDynamics=true;
		Model.ReadBinary(PopRef);
	}
	else
	{
		pCfgPtr->Sim.bRecycleDynamics=true;
		pCfgPtr->Sim.bDynamics=true;
	}

	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fDelta=pCfgPtr->ModelConstraints.Max.fBDecayDelta;

	int nIle(gnModelsPerParam);
	while(fDelta >= pCfgPtr->ModelConstraints.Min.fBDecayDelta )
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fBDecayDelta=fDelta;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		//cout<<Model.fTotalPenalty<<endl;
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.FullCycle();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fDelta<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fDelta/=1.05;
	}

}

//A sigma test for a bug
void MCMCInitialDistributionTest(const int nTestNum=1000000, const char *filename=NULL)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
	
	pCfgPtr->MCMCParameterUse.bLPowerLawGamma = true;
	pCfgPtr->MCMCParameterUse.bLPowerLawAlpha = true; 
	pCfgPtr->MCMCParameterUse.bLPowerLawBeta = true;
	pCfgPtr->MCMCParameterUse.bBDecayDelta = true;
	pCfgPtr->MCMCParameterUse.bBDecayKappa = true;

	pCfgPtr->MCMCParameterConstraint.bLPowerLawGamma = true;
	pCfgPtr->MCMCParameterConstraint.bLPowerLawAlpha = true; 
	pCfgPtr->MCMCParameterConstraint.bLPowerLawBeta = true;
	pCfgPtr->MCMCParameterConstraint.bBDecayDelta = true;
	pCfgPtr->MCMCParameterConstraint.bBDecayKappa = true;

	//read the config.xml
	if(filename != NULL) pCfgPtr->ReadXMLMCMC(filename);

	SModelParameters Params;
	CModelParameterDistributions Tester;
	const int nRes(1000);

	CHistogram Alpha(pCfgPtr->ModelConstraints.Min.fLPowerLawAlpha,
			 pCfgPtr->ModelConstraints.Max.fLPowerLawAlpha,
			 nRes);

	CHistogram Beta(pCfgPtr->ModelConstraints.Min.fLPowerLawBeta,
			 pCfgPtr->ModelConstraints.Max.fLPowerLawBeta,
			 nRes);

	CHistogram Gamma(log10(pCfgPtr->ModelConstraints.Min.fLPowerLawGamma),
			 log10(pCfgPtr->ModelConstraints.Max.fLPowerLawGamma),
			 nRes);

	CHistogram Delta(log10(pCfgPtr->ModelConstraints.Min.fBDecayDelta),
			 log10(pCfgPtr->ModelConstraints.Max.fBDecayDelta),
			 nRes);

	CHistogram Kappa(pCfgPtr->ModelConstraints.Min.fBDecayKappa,
			 pCfgPtr->ModelConstraints.Max.fBDecayKappa,
			 nRes);

	for(int i(0);i<nTestNum;i++)
	{
		Params = Tester.InitialDistribution();
		Alpha.AddPoint(Params.fLPowerLawAlpha);
		Beta.AddPoint(Params.fLPowerLawBeta);
		Gamma.AddPoint(log10(Params.fLPowerLawGamma));
		Delta.AddPoint(log10(Params.fBDecayDelta));
		Kappa.AddPoint(Params.fBDecayKappa);

		cout<<Params<<endl;
	}


	Alpha.Densify();
	Beta.Densify();
	Gamma.Densify();
	Delta.Densify();
	Kappa.Densify();



	Alpha.Norm();
	Beta.Norm();
	Gamma.Norm();
	Delta.Norm();
	Kappa.Norm();

	Alpha.Write("MCMCInitialDistributionTest.Alpha.hist");
	Beta.Write("MCMCInitialDistributionTest.Beta.hist");
	Gamma.Write("MCMCInitialDistributionTest.Gamma.hist");
	Delta.Write("MCMCInitialDistributionTest.Delta.hist");
	Kappa.Write("MCMCInitialDistributionTest.Kappa.hist");

}


void MCMCSamplingDistributionTest(const int nTestNum=1000000, const char *filename=NULL)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();

	pCfgPtr->MCMCParameterUse.bLPowerLawGamma = true;
	pCfgPtr->MCMCParameterUse.bLPowerLawAlpha = true; 
	pCfgPtr->MCMCParameterUse.bLPowerLawBeta = true;
	pCfgPtr->MCMCParameterUse.bBDecayDelta = true;
	pCfgPtr->MCMCParameterUse.bBDecayKappa = true;

	pCfgPtr->MCMCParameterConstraint.bLPowerLawGamma = true;
	pCfgPtr->MCMCParameterConstraint.bLPowerLawAlpha = true; 
	pCfgPtr->MCMCParameterConstraint.bLPowerLawBeta = true;
	pCfgPtr->MCMCParameterConstraint.bBDecayDelta = true;
	pCfgPtr->MCMCParameterConstraint.bBDecayKappa = true;

	//read the config.xml
	if(filename != NULL) pCfgPtr->ReadXMLMCMC(filename);

	SModelParameters CentreParams = pCfgPtr->Model;
	SModelParameters Params;
	CModelParameterDistributions Tester;

	const int nRes(1000);

	cout<<pCfgPtr->Model.fLPowerLawAlpha<<endl;
	cout<<pCfgPtr->Model.fLPowerLawAlpha - 3.*pCfgPtr->MCMCStepSize.fLPowerLawAlpha<<" "
	    <<pCfgPtr->Model.fLPowerLawAlpha + 3.*pCfgPtr->MCMCStepSize.fLPowerLawAlpha<<endl;

	CHistogram Alpha(pCfgPtr->Model.fLPowerLawAlpha - 3.*pCfgPtr->MCMCStepSize.fLPowerLawAlpha,
			 pCfgPtr->Model.fLPowerLawAlpha + 3.*pCfgPtr->MCMCStepSize.fLPowerLawAlpha,
			 nRes);
	CHistogram Beta(pCfgPtr->Model.fLPowerLawBeta -3.*pCfgPtr->MCMCStepSize.fLPowerLawBeta,
			 pCfgPtr->Model.fLPowerLawBeta +3.*pCfgPtr->MCMCStepSize.fLPowerLawBeta,
			 nRes);
	CHistogram Gamma(log10(pCfgPtr->Model.fLPowerLawGamma) -3.*log10(pCfgPtr->MCMCStepSize.fLPowerLawGamma),
			 log10(pCfgPtr->Model.fLPowerLawGamma) +3.*log10(pCfgPtr->MCMCStepSize.fLPowerLawGamma),
			 nRes);

	CHistogram Delta(log10(pCfgPtr->Model.fBDecayDelta) -3.*log10(pCfgPtr->MCMCStepSize.fBDecayDelta),
			 log10(pCfgPtr->Model.fBDecayDelta) +3.*log10(pCfgPtr->MCMCStepSize.fBDecayDelta),
			 nRes);
	CHistogram Kappa(pCfgPtr->Model.fBDecayKappa -3.*pCfgPtr->MCMCStepSize.fBDecayKappa,
			 pCfgPtr->Model.fBDecayKappa +3.*pCfgPtr->MCMCStepSize.fBDecayKappa,
			 nRes);

	for(int i(0);i<nTestNum;i++)
	{
		Params = Tester.SamplingDistribution(CentreParams);
		Alpha.AddPoint(Params.fLPowerLawAlpha);
		Beta.AddPoint(Params.fLPowerLawBeta);
		Gamma.AddPoint(log10(Params.fLPowerLawGamma));
		Delta.AddPoint(log10(Params.fBDecayDelta));
		Kappa.AddPoint(Params.fBDecayKappa);
		cout<<Params<<endl;
	}

	Alpha.Densify();
	Beta.Densify();
	Gamma.Densify();
	Delta.Densify();
	Kappa.Densify();


	Alpha.Norm();
	Beta.Norm();
	Gamma.Norm();
	Delta.Norm();
	Kappa.Norm();

	Alpha.Write("MCMCSamplingDistributionTest.Alpha.hist");
	Beta.Write("MCMCSamplingDistributionTest.Beta.hist");
	Gamma.Write("MCMCSamplingDistributionTest.Gamma.hist");
	Delta.Write("MCMCSamplingDistributionTest.Delta.hist");
	Kappa.Write("MCMCSamplingDistributionTest.Kappa.hist");

}










void TestA()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();
#if 0
//	pCfgPtr->ReadXML("GalacticPopulation.xml");
	pCfgPtr->pRan->m_PrepGaussCDF(0.01);
	//for(int i(0);i < (pCfgPtr->pRan->m_nGaussCDFSize) ;i++) cout<<pCfgPtr->pRan->m_pfGaussCDF[i]<<endl;
TTimestamp t0,t1,t2,t3;
	t0=GetTimestamp();
	for(int i(0);i < 10000000 ;i++) pCfgPtr->pRan->InverseTransformSamplingGauss();
	t1=GetTimestamp();
	for(int i(0);i < 10000000 ;i++) pCfgPtr->pRan->Gauss();
	t2=GetTimestamp();
	for(int i(0);i < 10000000 ;i++) pCfgPtr->pRan->GaussBruteForce(1.,0.);
	t3=GetTimestamp();
	cerr<<t1-t0<<endl;
	cerr<<t2-t1<<endl;
	cerr<<t3-t2<<endl;

for(int i(0);i < 10000000 ;i++) cout<<pCfgPtr->pRan->InverseTransformSamplingGauss()<<endl;
	return;

#endif
//	pCfgPtr->Sim.nMaxPulsarNum=gnMaxPulsarNum0;
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=false;
	pCfgPtr->Sim.fMaxTimeStep=0.1;
//	TestRandomizer();
	CModel Model;
	Model.ReadBinary("/home/mcie/Data/pop.bin");
	//CMCMC MCMC;
	RealType fLnL,fAvgLnL,fMinLnL,fMaxLnL;
	RealType fDelta=pCfgPtr->ModelConstraints.Max.fBDecayDelta;
	stringstream SS;
	string Pre;

	int nIle(gnModelsPerParam);
	while(fDelta >= pCfgPtr->ModelConstraints.Min.fBDecayDelta )
	{
		fAvgLnL=0.;

		pCfgPtr->Model.fBDecayDelta=fDelta;
		Model.Zero();
		Model.Evolve();
		Model.Densify();
		fLnL=Model.LnLikelihood();

		fAvgLnL+=fLnL;
		fMinLnL=fMaxLnL=fLnL;
		for(int i(0);i<nIle-1;i++)
		{	
			Model.Zero();
			Model.Evolve();
			Model.Densify();
			fLnL=Model.LnLikelihood();
			fAvgLnL+=fLnL;
			if(fMinLnL>fLnL) fMinLnL=fLnL;
			if(fMaxLnL<fLnL) fMaxLnL=fLnL;
		}
		fAvgLnL/=nIle;
		cout<<fDelta<<" "<<fAvgLnL<<" "<<fMinLnL<<" "<<fMaxLnL<<endl;
		fDelta/=1.05;
	}

}

void PrintHelp()
{
	cout<<"Model Tester"<<endl;
	cout<<"================================================================="<<endl;
	cout<<"Available options:"<<endl;
	cout<<"I.  Sequentional tests:"<<endl;
	cout<<"-SeqParam        - groud work for the sequentional test (sigmaBinit for the moment)"<<endl;
	cout<<"-A               - sequentional PowerLaw Alpha test"<<endl;
	cout<<"-B               - sequentional PowerLaw Beta test"<<endl;
	cout<<"-G               - sequentional PowerLaw Gamma test"<<endl;
	cout<<"-K               - sequentional Kappa test with NarayanL"<<endl;
	cout<<"-D               - sequentional Delta test with NarayanL"<<endl;
	cout<<"-PM              - sequentional PInitMean test with NarayanL"<<endl;
	cout<<"-PS              - sequentional PInitSigma test with NarayanL"<<endl;
	cout<<"-BM              - sequentional BInitMean test with NarayanL"<<endl;
	cout<<"-BS              - sequentional BInitSigma test with NarayanL"<<endl;
	cout<<"-Psi             - sequentional DeathPsi test with NarayanL"<<endl;
	cout<<"II.  MCMCTests:"<<endl;
	cout<<"-I               - initial parameters test"<<endl;
	cout<<"-S               - sampling test around the best parameters"<<endl;
	cout<<"III. CatalogueTests:"<<endl;
	cout<<"-C               - computing max dist/dm from dm/dist"<<endl;
	cout<<"IV.  Misc:"<<endl;
	cout<<"-h               - prints this help"<<endl;
	cout<<"-Dist            - PSR distributions"<<endl;
	cout<<"-ATNFCov         - ATNF Coverage"<<endl;
	cout<<"-Death           - Death Area test"<<endl;
	cout<<"-Erot            - Radio luminosity ~Erot test"<<endl;
	cout<<"-i <PopRef.bin>  - PopRef"<<endl;
//	cout<<"-o <Chains.hdf>  - MCMC-type output for sequentional tests (avg LnLikelihood in case o n>1)"<<endl;
//	cout<<"-p <PopOut.m.hdf> - output"<<endl;
//	cout<<"-x <Config.xml>  - Xml config"<<endl; 
	cout<<"-DM              - computing DMtab from NE2001 model"<<endl;
	cout<<"-TDMCont         - test time DMContainer Basic model"<<endl;
	cout<<"-TNEComp         - test time NE2001 comp Basic model"<<endl;
	cout<<"-Stab            - test stability"<<endl;
	cout<<"-DMDiff          - test difference bewteen NE2001 and DMContainer"<<endl;
	cout<<"Additional control (you may add them to the above selected flag):"<<endl;
	cout<<"-N               - Number of simulated pulsars (default 10k)"<<endl;
	cout<<"-n               - Number of models per parameters (default 1)"<<endl;
//	cout<<"-m               - if specified does \"m\" models from min to max in given parameter"<<endl;
	cout<<"================================================================="<<endl;
	cout<<"-Mki             - Hobbs' kicks (Maxwell 1D dist with a=265)"<<endl;
	cout<<"TODO:"<<endl;
//	cout<<"input xml"<<endl;
	cout<<"desc of output"<<endl;
}

int main(int argc, char ** argv)
{
	if(cmdOptionExists(argv, argv+argc, "-h") || argc==1) 
	{
		PrintHelp();
	}

	if(cmdOptionExists(argv, argv+argc, "-n"))
	{
		char * cNum = getCmdOption(argv, argv + argc, "-n");
		gnModelsPerParam=atoi(cNum);
	}

	if(cmdOptionExists(argv, argv+argc, "-N"))
	{
		char * cNum = getCmdOption(argv, argv + argc, "-N");
		gnMaxPulsarNum=atoi(cNum);
	}			

	char * input = getCmdOption(argv, argv + argc, "-i");


	if(cmdOptionExists(argv, argv+argc, "-A")) 
	{
		PowerLawLAlphaTest();
	}
	if(cmdOptionExists(argv, argv+argc, "-B")) 
	{
		PowerLawLBetaTest();
	}
	if(cmdOptionExists(argv, argv+argc, "-G")) 
	{
		PowerLawLGammaTest();
	}

	if(cmdOptionExists(argv, argv+argc, "-D")) 
	{
		NarayanLBDeltaSeqTest(input);
	}

	if(cmdOptionExists(argv, argv+argc, "-K")) 
	{
		NarayanLBKappaSeqTest();
	}
	if(cmdOptionExists(argv, argv+argc, "-R")) 
	{
		PowerLawLRotTest();
	}


	if(cmdOptionExists(argv, argv+argc, "-I"))
	{
		MCMCInitialDistributionTest();
	}
	if(cmdOptionExists(argv, argv+argc, "-S"))
	{
		MCMCSamplingDistributionTest();
	}
	if(cmdOptionExists(argv, argv+argc, "-C")) 
	{
		TestMaxDistATNF();
	}

	if(cmdOptionExists(argv, argv+argc, "-DM")) 
	{
		TestDMContainer();
	}
	
	if(cmdOptionExists(argv, argv+argc, "-PM")) 
	{
		NarayanLPInitMeanTest();
	}
	if(cmdOptionExists(argv, argv+argc, "-PS")) 
	{
		NarayanLPInitSigmaTest();
	}
	if(cmdOptionExists(argv, argv+argc, "-BM")) 
	{
		NarayanLBInitMeanTest();
	}
	if(cmdOptionExists(argv, argv+argc, "-BS")) 
	{
		NarayanLBInitSigmaTest();
	}

	if(cmdOptionExists(argv, argv+argc, "-Psi")) 
	{
		NarayanLDeathPsiTest();
	}


	if(cmdOptionExists(argv,argv+argc, "-TDMCont"))
	{
		TestModelTimeDMContainer();
	}

	if(cmdOptionExists(argv,argv+argc, "-TNEComp"))
	{
		TestModelTimeNE2001();
	}

	if(cmdOptionExists(argv,argv+argc, "-Stab"))
	{
		TestModelStability();
	}

	if(cmdOptionExists(argv,argv+argc, "-DMDiff"))
	{
		TestDMDiff();


	}
	
	if(cmdOptionExists(argv,argv+argc, "-Dist"))
	{
		TestDistributions();


	}

	if(cmdOptionExists(argv,argv+argc, "-ATNFCov"))
	{
		TestATNFCoverage();
	}


	if(cmdOptionExists(argv,argv+argc, "-Death"))
	{
		TestPulsarDeathArea();
	}
	if(cmdOptionExists(argv,argv+argc, "-Erot"))
	{
		TestLErot();
	}

	if(cmdOptionExists(argv,argv+argc, "-SeqParam"))
	{

		SeqentionaParameterTest();
	}

	if(cmdOptionExists(argv,argv+argc, "-Mki"))
	{

		TestMaxwell();
	}



	return 0;
}


