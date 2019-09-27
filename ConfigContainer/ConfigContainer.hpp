//FIXME poprawic naze plikow aby lepiej odpowiadala klasie

#pragma once
#include<Definitions.hpp>
#include<Randomizer.hpp>
#include<ModelParameters.hpp>
#include<Survey.hpp>
#include<Catalogue.hpp>
#include<CatalogueMasks.hpp>
#include<pugixml.hpp>
#include<DMContainer.hpp>

enum EKickingMethod
{
	PaczynskiKicking = 0,
	HobbsKicking = 1
//	ChuckNoris = 2
};



enum EExtremumMethod
{
	Likelihood=0,
	KolmogorovSmirnov=1
};

enum ELuminosityModel
{
	Narayan=0,
	Proszynski=1,
	PowerLaw=2,
	Random=3,
	Erot=4
};


enum EModelDensityContainer
{
	Density1x3D = 0,
	Density2x3D = 1,
	Density1x3DSparse = 2,
	Density2x3DSparse = 3,
	Density4DSparse = 4,
	Density3x1D = 5
};


struct SSimulationParameters 
{
	RealType fLnLikelihoodSigma;
	int nNumberOfIteration;///<Number of physical model iterating over the same geometrical model. A fix to boost the efficiency of the model.
	EKickingMethod KickingMethod;
	EExtremumMethod ExtremumMethod; ///<Extremum method 0 - Likelihood, 1 - Kolmogorov-Smirnov
	RealType fMaxPulsarAge; ///< [Myr] The maximal age (evolution time) for a pulsar given random init parameters
	RealType fMaxTimeStep;  ///< [Myr] The maximal timestep that produce good results (constrainted by Position Verlet method)

	RealType fMaxDist;      ///< [kpc] Maximal distance form the Earth-Observer.

	int nNumOfChains;	///< Number of chains in MCMC
	int nNumPerChain;	///< Number of links per chain in MCMC
	bool bSurveyGeometry;   ///< Flag if the survey should check for sky l,B 
	bool bCutATNFOutsideArea;

	bool bRecycleDynamics;  ///< Flag for reusing evolved spatial distribution (for speedup)
	bool bDynamicsComputed; ///< Flag to use in conjunction with bRecycleDynamics 
	bool bObservations;     ///< Flag for turning on the survey
	bool bDynamics;         ///< Flag for turning on the dynamic evolution insinde the galactic potential

	bool bPhysics;          ///< Flag for turning on the physical evolution of the pulsar (B, P, Pdot)

	bool bVerbose;          ///< Flag for the verboisty of the output
	bool bMCMCVerbose; 	///< Not as detailed as bVerbose by prints the MCMC sim output
	bool bOnlyPowerLaw; ///<Flag for MCMC: disable random, only power law Luminoisty

  
	int nMaxPulsarNum;      ///< The maximal number of pulsars //FIXME czy to jest uzywane??
	int nMaxAvailablePulsarNum; ///<The maximal number of pulsars that can be processed
	bool bCheck4Relevance;  ///<Enable checking if a pulsar can be seen from any of the surveys
	bool bReuseDynamics;	///< If Likelihood computation should reuse evolved spatial distribution (due to the DM computation cost)
	SDataFields DataFields; ///< Maska na obecnosc danych w katalogu
	RealType fDensityMinValue;

	//RealType fIntervalFraction;///<Parameters' jump interval of [min,max] in MCMC
	//RealType fIntervalScalingFactor;///<Scaling of the jump inteval trough tha MCMC chain
	RealType fGaussianDensitySigma;
	bool bGaussianDensityMod;
	EModelDensityContainer ModelDensityContainer;///<Choice of the method/container for the model's probability density: 1 - Density2x3D, 3 - Density2x3DSparse, 4 - Density4DSparse, 5 - Density3x1D 
	ELuminosityModel LuminosityModel;///<Choice of Luminosity model: 0 - Narayan, Ostriker, 1 - Proszynski, Przybycien, 2 - PowerLaw, 3 - Random, 4 - Erotational

	char sCatalogueFile[256];
	bool bUseCatalogue;///<A flag for acticating ATNF in postConfigInit
	bool bStrictLikelihood;///<If a model doesn't cover even one obs it has min value
	bool bATNFDeathLinesCut;///<Cuts out the cotalouge pulsars that fall under the death lines

	char sDMDataFile[256];
	bool bUsePrecomputedDMData;
	bool bUseDMContainer;
	long int nDMDistSlices;
	long int nHEALPixNSide;
	//RealType fMaxDist;
	bool bComputeDMDiff;///<compute difference betwean tabbed dm and NE2001 instead of DM.

	bool bComputeCatalogueMaxDist;///<A flag to compute ATNF max dist fMaxDist used if false

	bool bVerboseParameterUsage;///<Flag mainly for debugging parameters in MCMC
	bool bMCMCSequential;///>Do seqentional 1D chains instead of MCMC.

	bool bLnLikelihoodPenalty;///>Use penalization scheme
	RealType fLnLikelihoodPenaltyCoefficient;///>The weight of the penalty


	bool bPhysicalDensity;///>If the Density container has dV in physical units (or log of them)


	bool bSameSeedInMCMCModels;///>Will cause to reseed each time a model is started. Otherwise the prvided seed (as option or in xml) is used at the start of the computation (node-wise in case of cluster).

	bool bLEfficiency;///>Radio Luminosity efficiency as described at Szary at al.
	bool bLogLCorr;///>Random correction to the radio luminosity after Bates and Faucher-Giguere


	bool bNoDeathLine;

	RealType fDetectionTreshold;///>The orignal value was 10. - from inital consulatation with Stefan.
    bool bScaleJumps;
    RealType fJumpScallingFactor;
    int nScaleEveryJumps;

    bool bStartMCMCFromModel;
    bool bContinueMCMCFromLastPostion;
    int nLastMCMCPostion;
    int nLastMCMCChain;

    bool bUseRandomParamSampler;
};

struct SOutputParameters 
{
	int nDummy;

};


struct SPlotParameters
{
	///P-Pdot-L
	int nResLogP;
	int nResLogPdot;
	int nResLogS1400;

	RealType fZeroLogP;
	RealType fZeroLogPdot;
	RealType fZeroLogS1400;

	RealType fMaxLogP;
	RealType fMaxLogPdot;
	RealType fMaxLogS1400;


	RealType fMaxDM;
	int nResGalL;
	int nResGalKsi;
	int nResDM;


};
//FIXME to potrzebuje jakiegos prostszego rozwiazania - niby kto bedzie odpowiedzialny za sprzatanie pID?
struct SSurveySelection
{
	ESurveyID *pID;
	int	nNum;
	SSurveySelection()
	{
		nNum=0;
	};
	~SSurveySelection()
	{
		if(nNum) delete [] pID;
	};


};

struct SSurveysCuts
{
	bool bInitiated;
	RealType *pMaxSurveyDist;
	SSurveysCuts()
	{
		bInitiated=false;
	}
};

//Klasyczny leniwy singleton
class CConfiguration
{
public:
	

	SSimulationParameters 	Sim;
	SModelParameters 	Model;
	SModelParameters 	MCMCStepSize;
	SModelConstraints	ModelConstraints;
	SPlotParameters		Plot;
	SOutputParameters 	Out;
	SSurveySelection	SurveySelection;
	SSurveysCuts            SurveysCuts;

	SModelParametersMask 	MCMCParameterUse;
	SModelParametersMask 	MCMCParameterConstraint;

	
	CRandomizer 	*pRan;
	CCatalogue	*pATNF;
	CDMContainer    *pDM;
	

	static CConfiguration* GetInstance();
	void Print();
	void Reseed(const unsigned int nSeed, const unsigned nCalls=0);

	void PostConfigInit();
	//void BasicConfig();
	void ReadXMLMCMC(const char *FileName,const bool bVerbose=false);
	void ReadXML(const char *FileName);
	void WriteXMLMCMC(const char *FileName);
	void WriteXML(const char *FileName);
	~CConfiguration();

private:
	void BasicConfig();
	bool bCreatedDM;
	CConfiguration();                   
        void AppendMCMC(pugi::xml_node *pConfig);
        void AppendModel(pugi::xml_node *pConfig);
	void ReadMCMC(pugi::xml_document *pDoc, const bool bVerbose=false);
	void ReadModel(pugi::xml_document *pDoc, const bool bVerbose=false);

	// Dont forget to declare these two. You want to make sure they
        // are unaccessable otherwise you may accidently get copies of
        // your singleton appearing.
        CConfiguration(CConfiguration const&);              // Don't Implement
        void operator=(CConfiguration const&); // Don't implement
};

