#pragma once


#include<Randomizer.hpp>
#include<Definitions.hpp>
#include<Coordinates.hpp>
#include<Statistics.hpp>
#include<BasicStar.hpp>


RealType MYears2Seconds(const RealType fMYears);
RealType Seconds2MYears(const RealType fSeconds);



class CPulsarPhysics : public virtual CBasicStar
{
protected:
	CRandomizer *m_pRanPtr;
	CConfiguration *m_pCfgPtr;

	RealType m_fBFiledTimeDelatDistribution();
	RealType m_fBFieldTimeDeltaGaussDistribution(const RealType fMean, const RealType fSigma);

//	RealType m_fBFieldTimeDelta;
	RealType m_FinalLogBFieldDistribution();
	RealType m_InitLogBFieldDistribution();
	RealType m_InitPeriodDistribution();
	RealType m_InitLogPeriodDistribution();

	RealType m_PeriodDerivative(const RealType fP, const RealType fB);
	RealType m_NextBField(const RealType fB, const RealType fTimeStep);

	void     m_EvolveBField(const RealType fTimeStep);
	RealType m_RombergIntegralBFieldSquared(const RealType fTimeStep);

	void m_ComputePPdotBField(const RealType fTimeStep);
	RealType m_PeriodFunction(const RealType fTimeStep);

public:

	//RealType fFinalBField; //Przeniesione do CBasicStar
	CPulsarPhysics();
	~CPulsarPhysics();
	void Init();
	void Evolve(const RealType fTimeStep, const bool bOneStep=false);
	void EvolveRK4(const RealType fTimeStep);
	void EvolveRK4Adaptive(const RealType fMaxTimeStep);
	void InitPulsarPhysics();
	void EvolveForwardNewton(const RealType fTimeStep);
};


class CPulsarPhysicTester : public CPulsarPhysics
{
public:
	void InitPeriod(const int nNum, const char* FileName);
	void FinalBField(const int nNum, const char* FileName);
	void InitBField(const int nNum, const char* FileName);
};


