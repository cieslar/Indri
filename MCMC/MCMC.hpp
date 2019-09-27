#pragma once
#include<Definitions.hpp>

#include<Randomizer.hpp>
#include<Statistics.hpp>

#include<fstream>
#include<iostream>

#include<cmath>

using namespace std;

#include<ModelParameters.hpp>
#include<ConfigContainer.hpp>

#include<Population.hpp>





///Class that glue togheter observations and synthetic population
class CModel
{
protected:
	CConfiguration	*m_pCfgPtr;

	CPopulation	*m_pPop;///<This holds the exact pulsars' data.
	CObservations	**m_ppObs;///<Subsets of the catalogue (for each Survey)
	CDensity	**m_ppDensity;///<Model's density coresponding to each survey
	CRadioDetector  *m_pRadioDet;///<A patch to access the visibility functions.
	int 		m_nObsNum;///<Number of considered surveys
	

	RealType m_fWorstCaseLnLikelihood;


public:
	void Densify();
	RealType LnLikelihood();
	RealType DValue();

	void ComputeMinimizationValue();///<Computes both LnLikelihood and DValue.

	RealType fTotalPenalty;

	RealType Efficiency;

	void Evolve();///Evolves the model form the start to end with parameters form pCfgPtr
	void Zero();///Reset everything
	//void ZeroPhysics();///Resets only the Physics part of the evolved PSR. The dynamics (position in Galaxy is unchanged)
		
	CModel();
	//CModel(const char *FileName);///Constructor with reading existing population
	~CModel();
	//Przydadzą się też funkcje wypisujące przekroje, itd.
	//Ale poki co wersja minimalna tak aby dalo sie ja podpiac pod
	void Write(const char* Prefix=NULL);
	void ReadBinary(const char* PopFile);
	void ReadBinary(const char* PopFile, CConfiguration *pCfgPtr);

	void WriteBinaryHDF5(const char* FileName,const int nIteration);
	void WritePartialBinaryHDF5(const char* FileName);
	void WriteAsCatalogue(const char* FileName, const int nNumber);

	void PrintPulsar(const int nNum) const;
	void Test();
	void Check();


	void FillDensityContainers();
	void NormalizeDensityContainers();
	void ZeroPop();
	void FullCycle(const char* HDFFileName = NULL);
};


class CTestModel: public CModel
{
public:
	CTestModel(){};
	~CTestModel(){};

	void Test();
};


#if 0


//jak jest ogniwo??
struct Link
{
	SModelParameters Point;
	RealType fValue;
};

RealType Fun(const SModelParameters P)
{
	//return P.fX * P.fY * exp(-P.fX) * exp(-P.fY) * sin(2.*M_PI/P.fX) * cos(P.fY);
	
	//if(P.fX < 0. || P.fX >10.) return 0.;
	//return P.fX*exp(-P.fX);
	return 0.;
	
}
//jak bede wyliczal likelihood - jak uwzglednic blad pomiarowy? rozklad normalny? rozklad jednorodny?
//czy da sie zrobic 

#endif
class CModelParameterDistributions
{
private:
	CConfiguration *m_pCfgPtr;
	CRandomizer *m_pRanPtr;
	RealType m_fIntervalFraction;

public:
	CModelParameterDistributions()
	{
		m_pCfgPtr = CConfiguration::GetInstance();
		m_pRanPtr = m_pCfgPtr->pRan;
		m_fIntervalFraction = m_pCfgPtr->ModelConstraints.fIntervalFraction;

	};
	~CModelParameterDistributions(){};

	SModelParameters InitialDistribution() const; ///< Probability distribution for starting the chain

	SModelParameters SamplingDistribution(const SModelParameters Current) const; ///< Probability distribution for making a link in the chain
	SModelParameters BasicSamplingDistribution(const SModelParameters Current) const; ///< Probability distribution for making a link in the chain
	SModelParameters RandomParamSamplingDistribution(const SModelParameters Current) const; ///< Probability distribution for making a link in the chain

	void ScaleJumpSigmas();
protected:

	//PowerLaw
	RealType m_PowerLawGammaInitialDistribution() const; 
	RealType m_PowerLawAlphaInitialDistribution() const; 
	RealType m_PowerLawBetaInitialDistribution()  const; 

	RealType m_PowerLawGammaSamplingDistribution(const RealType CurrentGamma) const;
	RealType m_PowerLawAlphaSamplingDistribution(const RealType CurrentAlpha) const;
	RealType m_PowerLawBetaSamplingDistribution (const RealType CurrentBeta)  const;


	//BDecay
	RealType m_BDecayDeltaInitialDistribution() const;
	RealType m_BDecayKappaInitialDistribution() const;

	RealType m_BDecayDeltaSamplingDistribution(const RealType CurrentDelta) const;
	RealType m_BDecayKappaSamplingDistribution(const RealType CurrentKappa) const;

	RealType m_BInitMeanInitialDistribution() const;
	RealType m_BInitSigmaInitialDistribution() const;

	RealType m_BInitMeanSamplingDistribution(const RealType CurrentMean) const;
	RealType m_BInitSigmaSamplingDistribution(const RealType CurrentSigma) const;

	RealType m_PInitMeanInitialDistribution() const;
	RealType m_PInitSigmaInitialDistribution() const;

	RealType m_PInitMeanSamplingDistribution(const RealType CurrentMean) const;
	RealType m_PInitSigmaSamplingDistribution(const RealType CurrentSigma) const;

	RealType m_DeathPsiInitialDistribution() const;
	RealType m_DeathPsiSamplingDistribution(const RealType CurrentPsi) const;



};

class CMCMC: protected CModelParameterDistributions
{
protected:
	int m_nNumOfChains;
	int m_nNumPerChain;
	bool m_bChainMemOn;


	CConfiguration *m_pCfgPtr;
	CRandomizer *m_pRanPtr;

	SModelParameters **m_pChain; ///< Container for the Markov chains

	SModelParameters m_Params; ///< Container for the current model's parameters

	CModel *m_pModel; ///< The main computational workhorse 

	void m_CreateEmptyChains();
	void m_DeleteChains();

	void m_ComputeValue(SModelParameters& Parameters);
	SModelParameters m_StartChain();
	SModelParameters m_MakeLink(SModelParameters Previous);
	bool m_CheckParameterOutsideBounds(const SModelParameters Params) const; ///< returns true if parameters are outside predefined bounds.
    void m_VisitLastPositionInChain(const int nChain, const int nLastLink);
public:

	void Test();
	void MakeChains();
	void Write();
	void Write(const int nPointOfChain);
	void WriteWholeChains();
	void WriteBinary(const char* Filename);
	void WriteBinaryHDF5(const char* Filename) const;
	void CombineBin(const char* OutputFilename, char** InputFilenameVec, const int nNumOfFiles);
	void ReadBinary(const char* Filename);
	void PopReadBinary(const char* PopFile);
	void MakeDensitySpace();

	CMCMC();

	~CMCMC();

};

class CMCMCPostProcessing : public CMCMC
{
	SModelParameters MinValues;
	SModelParameters MaxValues;

public:
	CMCMCPostProcessing(){};
	~CMCMCPostProcessing(){};
	void SetExtremes();
	void FindExtremes();
	void WriteExtremes();
	void MakeHistograms(const int nSquareResolution=100);
	void MakeHistograms1D(const int nResolution=100);
};



