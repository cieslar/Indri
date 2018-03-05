#pragma once

#include<ConfigContainer.hpp>
#include<PulsarDynamics.hpp>
#include<PulsarPhysics.hpp>

#include<Projections.hpp>
#include<Definitions.hpp>

#include<NE2001.hpp>
#include<string>

class CPulsar: public CPulsarDynamics, public CPulsarPhysics
{
protected:
	bool m_CheckDistance();
	bool m_CheckBeaming();
	bool m_Check4Relevance() const;///< checks if the pulsar can be seen by any of the surverys;
	CConfiguration *m_pCfgPtr;
	CRandomizer *m_pRanPtr;
	SBasicStar *m_pBS;
public:

	CPulsar();
	CPulsar(const SBasicStar& Star);
	~CPulsar();
	void Evolve(const RealType fMaxTimeStep);
	void EvolveDynamics(const RealType fMaxTimeStep);
	void EvolvePhysics(const RealType fMaxTimeStep);
	RealType LogRadioLuminosity();
	RealType LogL400();

	RealType LogL1400Random() const;///<[mJy*kpc^2]
	RealType LogL400PowerLaw() const;///<[mJy*kpc^2]

	RealType SpectralIndexDistribution() const;
	bool CheckProbability() const;

	RealType DeathLinesProbability() const;
	RealType LogL400Erot() const;

	RealType L400Proszynski() const;///<[mJy*kpc^2]
	RealType RadioFlux400() const; ///result in [mJy], freq in [MHz]
	RealType RadioFlux(const RealType fFrequency) const; ///result in [mJy], freq in [MHz]
	
	RealType BeamingFactor(const RealType fP) const; ///<as a probability [0:1]
	RealType BeamingFactor() const; ///<as a probability [0:1]

	RealType RadioLuminosity400(const RealType fLuminosity,const RealType fFrequency ) const;

	bool AboveDeathLines(const RealType fP, const RealType fPDot) const;
	bool AboveDeathLines() const;

	RealType  LogPFromDeathLine(const RealType fLogPDot) const;
	RealType  LifeProbability(const RealType fLogP, const RealType fLogPDot, const RealType fPhi) const;

	void Copy(const CPulsar *PSR);
	RealType DispersionMeasure();
	void Project();
	void Init(const RealType fPositionStartTime);
	//void Print() const;

	CPulsar& operator=(const SBasicStar& BS);

	RealType LEfficiency() const;
	RealType LogLCorr() const;

	std::string AsCatalogueEntry(const int nRefNumber) const;
	RealType ETot();
};

class CPulsarTester: public CPulsar
{
public:
	void  DeathTester() const;
	void  LErotTester(const RealType fGamma=pow(10.,1.635), const RealType fAlpha=1.);
    void CatModelLTest() ;
};
