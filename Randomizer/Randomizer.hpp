#pragma once

#include<Definitions.hpp>
#include<Coordinates.hpp>
#include<MersenneTwister.h>
#include<Debug.hpp>

enum RandomizerTypeEnum {DefaultRandomizer=0, MTRandomizer=1};

void TestRandomizer();

class CRandomizer
{
public:
	unsigned int m_nSeed;
	unsigned long int m_nCalls;
	unsigned long int m_nCallsInt;//nie sprawdzilem czy jest osobny mechanizm losowania dla intow stad osobne
	RealType InverseTransformSamplingGauss(const RealType fSigma=1., const RealType fMi=0.);

	RealType m_RombergIntegralGauss(const RealType fFrom, const RealType fTo);
	void m_PrepGaussCDF(const RealType fStep);

	int m_nGaussCDFSize;
	bool m_bGaussCDF;
	RealType m_fGaussCDFStep;
	RealType *m_pfGaussCDF;


	const int m_nSigmaRange=6;
	const RealType m_fSigmaToleranceInterval=0.999999998027;
	const RealType m_fGaussCDFResolution=0.001;


	RealType m_MaxwellDistribution(const RealType fX, const RealType fA);
	RealType Maxwell(const RealType fA);

	MTRand* m_pMTR;
	RealType m_GaussDistribution(const RealType fX, const RealType fSigma=1., const RealType fMi=0.);
	RealType m_ExponentialDistribution(const RealType fX, const RealType fLambda);
	RealType m_BimodalGaussDistribution(const RealType fX, const RealType fWeight, const RealType fSigma1, const RealType fSigma2, const RealType fMi1=0., const RealType fMi2=0.);
//	RealType m_LogNormalDistribution(const RealType fX, const RealType fSigma, const RealType fMi=0.);
	void Seed(const unsigned int nSeed, const unsigned int nRollCount=0 );
//public:
	CRandomizer();
	~CRandomizer();

	RealType Random();
//	int RandomInt();
	bool RandomLogical();
	RealType Flat();///< Simple flat distribution [0,1.]
	RealType FlatOpen();///< (0,1)
	RealType FlatLeftOpen();///< (0,1]
	RealType FlatMinMax(RealType fMin, RealType fMax);
	RealType LogFlatMinMax(RealType fMin, RealType fMax);
	RealType FlatMeanSigma(RealType fMean, RealType fSigma);
	RealType LogFlatMeanSigma(RealType fMean, RealType fSigma);
	RealType PositiveFlatMeanSigma(RealType fMean, RealType fSigma);
	RealType PositiveLogFlatMeanSigma(RealType fMean, RealType fSigma);

	RealType LogFlatMeanLogSigma(RealType fMean, RealType fLogSigma);
	RealType PositiveLogFlatMeanLogSigma(RealType fMean, RealType fLogSigma);

	RealType GaussBruteForce(const RealType fSigma, const RealType fMi=0.);
	RealType Gauss(const RealType fSigma=1., const RealType fMi=0.);
	RealType GenerateGaussianBoxMuller(const RealType fSigma, const RealType fMi=0.);
	Cartesian Gauss2D(const RealType fSigma, const RealType fMi=0.);
	RealType BimodalGauss(const RealType fWeight, const RealType fSigma1, const RealType fSigma2, const RealType fMi1=0., const RealType fMi2=0.);

	RealType PositiveGauss(const RealType fSigma, const RealType fMi=0.);
	RealType PositiveGaussBruteForce(const RealType fSigma, const RealType fMi=0.);	

	RealType Exponential(const RealType fLambda);
	RealType DoubleSidedExponential(const RealType fLambda);
//	RealType LogNormal(const RealType fSigma, const RealType fMi=0.);

	void PrintInfo();
	
	TTimestamp nRanCount, nMathCount;
	TTimestamp nTotCalls;

	Cartesian UniformShpereDirection();
};

void TestRandomizer();
