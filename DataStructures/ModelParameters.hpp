#pragma once

#include<Definitions.hpp>
#include<map>
#include<iostream>
#include<string>
#include<CatalogueMasks.hpp>

using namespace std;

struct SModelParameters
{
	unsigned int nSeed;


	RealType fLPowerLawGamma;
	RealType fLPowerLawAlpha;
	RealType fLPowerLawBeta;
	RealType fBDecayDelta;
	RealType fBDecayKappa;//wykladnik w exp(x^kappa/Delta)

	RealType fBInitMean;///<log10(B/1[GS])
	RealType fBInitSigma;///<log10(B/1[GS])
	RealType fPInitMean;///<log10(P/1[s])
	RealType fPInitSigma;///<log10(P/1[s])
	RealType fVKickMod1Mean;///<just a place holder not implemented
	RealType fVKickMod1Sigma;///<just a place holder not implemented
	RealType fVKickMod2Mean;///<just a place holder not implemented
	RealType fVKickMod2Sigma;///<just a place holder not implemented

	RealType fBFinalMean;///<log10(B/1[GS])
	RealType fBFinalSigma;///<log10(B/1[GS])

	

	RealType fDeathPsi;///<Parameter for the Death Area (width of the atan distr.)
//	RealType fBDecayTimeDeltaMean;//[Myr]
//	RealType fBDecayTimeDeltaSigma;//[Myr]

	RealType fLnLikelihood;
	RealType fDValue;
	//bool bLPowerLaw;
	//bool bLRandom;



 //dlaczego nie mogę mieć prostego konstruktora u użyć tego potem w unii?
	SModelParameters()
	{
//		bLPowerLaw=false;
//		bLRandom=false;
	//	fBDecayTimeDeltaMean=0.;
	//	fBDecayTimeDeltaSigma=0.;
		fLnLikelihood=0.;
		fDValue=0.;
		nSeed=0;


		fLPowerLawGamma=0.;
		fLPowerLawAlpha=0.;
		fLPowerLawBeta=0.;

		fBDecayDelta=0.;
		fBDecayKappa=0.;//wykladnik w exp(x^kappa/Delta)

		fBInitMean=0.;	
		fBInitSigma=0.;
		fBFinalMean=0.;	
		fBFinalSigma=0.;	
	
		fPInitMean=0.;
		fPInitSigma=0.;
		fVKickMod1Mean=0.;
		fVKickMod1Sigma=0.;
		fVKickMod2Mean=0.;
		fVKickMod2Sigma=0.;
		fDeathPsi=0.;
	};


	//FIXME znalezc sposob na zerowanie
};


std::ostream& operator<<(std::ostream& out, const SModelParameters& Model);
std::string DescribeModelParameters();

struct SModelParametersMask
{
	bool bLPowerLawGamma;
	bool bLPowerLawAlpha;
	bool bLPowerLawBeta;
	bool bBDecayDelta;
	bool bBDecayKappa;//wykladnik w exp(x^kappa/Delta)

	bool bBInitMean;	
	bool bBInitSigma;
	bool bBFinalMean;	
	bool bBFinalSigma;		
	bool bPInitMean;
	bool bPInitSigma;
	bool bVKickMod1Mean;
	bool bVKickMod1Sigma;
	bool bVKickMod2Mean;
	bool bVKickMod2Sigma;
	bool bDeathPsi;


	SModelParametersMask()
	{
		bLPowerLawGamma=false;
		bLPowerLawAlpha=false;
		bLPowerLawBeta=false;
		bBDecayDelta=false;
		bBDecayKappa=false;
		bBFinalMean=false;	
		bBFinalSigma=false;	

		bBInitMean=false;	
		bBInitSigma=false;	
		bPInitMean=false;
		bPInitSigma=false;
		bVKickMod1Mean=false;
		bVKickMod1Sigma=false;
		bVKickMod2Mean=false;
		bVKickMod2Sigma=false;
		bDeathPsi=false;
	};
};



#if 0
#define MODEL_PARAMETERS_SIZE 24

union UModelParameters
{
	SModelParameters Data;
	char	Word[MODEL_PARAMETERS_SIZE];

};

#endif
//Opisuje zakres losowania rozkladow pierwotnych
//pierwsze przyblizenie to rozklad plaski o brzegach Min i Max
struct SModelConstraints
{
	///Flags to activate the variables as a model's parameters
	bool bBDecayTimeDeltaMean;
	bool bBDecayTimeDeltaSigma;

	SModelParameters Min;
	SModelParameters Max;	
	//SModelParameters Sigma;	
	RealType fIntervalFraction;///<What fraction of [Min,Max] should be the MCMC jump
};

struct SProbabilitySpacePoint4D
{
	int nLogP;
	int nLogPdot;
	int nLogS1400;
	int nDM;
	//int nHEALPix;
};

struct SCenterOfMass4D
{
	RealType fLogP;
	RealType fLogPdot;
	RealType fLogS1400;
	RealType fLogDM;
	//int nHEALPix;
};




bool operator<(const SProbabilitySpacePoint4D &A, const SProbabilitySpacePoint4D &B);

typedef map<SProbabilitySpacePoint4D, RealType>::iterator IteratorMapType;



//FIXME przepisac na "A"
struct SProbabilitySpacePointRadio
{
	int nLogP;
	int nLogPdot;
	int nLogS1400;
};

struct SProbabilitySpacePointRadioSmooth
{
	RealType fLogP;
	RealType fLogPdot;
	RealType fLogS1400;
};




struct SCenterOfMassRadio
{
	RealType fLogP;
	RealType fLogPdot;
	RealType fLogS1400;
};


struct SProbabilitySpacePointSpace
{
	int nGalL;
	int nGalKsi;
	int nDM;
};






bool operator<(const SProbabilitySpacePointRadio &A, const SProbabilitySpacePointRadio &B);
bool operator<(const SProbabilitySpacePointSpace &A, const SProbabilitySpacePointSpace &B);

typedef map<SProbabilitySpacePointRadio, RealType>::iterator IteratorMapTypeRadio;
typedef map<SProbabilitySpacePointSpace, RealType>::iterator IteratorMapTypeSpace;

SModelParameters operator+(const SModelParameters &A, const SModelParameters &B);
SModelParameters operator*(const SModelParameters &A, const RealType &fB);
SModelParameters operator*(const RealType &fB, const SModelParameters &A);
SModelParameters operator-(const SModelParameters &A, const SModelParameters &B);
SModelParameters pow(const SModelParameters A, const RealType fB);
SModelParameters log10(const SModelParameters A);
SModelParameters sqrt(const SModelParameters A);


SProbabilitySpacePointSpace ComputeSpacePointCoordiantesSpace(const SCatalogueEntry PSR);
SProbabilitySpacePointRadio ComputeSpacePointCoordiantesRadio(const SCatalogueEntry PSR);
SProbabilitySpacePointSpace ComputeSpacePointCoordiantesSpace(const CPulsar PSR);
SProbabilitySpacePointRadio ComputeSpacePointCoordiantesRadio(const CPulsar PSR);

SProbabilitySpacePointRadioSmooth ComputeSpacePointCoordiantesRadioSmooth(const SCatalogueEntry PSR);
SProbabilitySpacePointRadioSmooth ComputeSpacePointCoordiantesRadioSmooth(const CPulsar *pPSR);

RealType GaussianishDensityContribution(const SProbabilitySpacePointRadio RefSpacePoint, const SProbabilitySpacePointRadioSmooth SpacePoint, const RealType fSigma);



SProbabilitySpacePoint4D ComputeSpacePointCoordiantes4D(const CPulsar PSR);
SProbabilitySpacePoint4D ComputeSpacePointCoordiantes4D(const SCatalogueEntry PSR);
