#pragma once
#include<Definitions.hpp>
#include<Coordinates.hpp>
#include<ostream>
#include<H5Cpp.h>

struct SBasicStar
{
	bool bBeaming;///< flag for beaming true/false 
	bool bInsideGalaxy;///< flag for beeing inside Galaxy
	bool bIsLost;//nie brana dalej pod uwage	

	RealType fProbability;//Probability of detection due to different causes

	RealType fAge;//[Myr]
	RealType fRelativeAge;///< Relative to observation timepoint.[Myr]
	RealType fMass;//[MSun]

	RealType fPeriod;//okres obroty
	RealType fPeriodDot;//pochodna okresu obrotu
	RealType fBField;//Indukcja pola magnetycznego
	RealType fInitialPeriod;//okres obroty
	RealType fInitialPeriodDot;//pochodna okresu obrotu
	RealType fInitialBField;//Początkowe pole mag
	RealType fFinalBField;///<Koncowe pole mag
	RealType fSpectralIndex;

	RealType fL400;///< [mJy*kpc^2] @400Mhz
	RealType fS400;///< [mJy] @1400 MHz 
	RealType fS1400;

	RealType fDM;///< DispersionMeasure

	Cartesian Position;//FIXME usunac m_ ///<(x,y,z)[kpc] 
	Cartesian Velocity;///<(Vx,Vy,Vz)[kpc/Myr]
	Cartesian InitialPosition;//FIXME usunac m_ ///<(x,y,z)[kpc] 
	Cartesian InitialKick;///<(Vx,Vy,Vz)[kpc/Myr]
	Galactic GalacticPos;///<(l,b,D) [rad],[rad],[kpc]
	Galactic GalacticPosDeg;///<(l,b,D) [deg], [deg], [kpc]
	

};

struct SPopulationDataSet
{
	H5::DataSet bBeaming;
	H5::DataSet bInsideGalaxy;
	H5::DataSet bIsLost;
	H5::DataSet fProbability;
	H5::DataSet fAge;
	H5::DataSet fRelativeAge;
	H5::DataSet fMass;
	H5::DataSet fPeriod;
	H5::DataSet fPeriodDot;
	H5::DataSet fBField;
	H5::DataSet fInitialPeriod;
	H5::DataSet fInitialPeriodDot;
	H5::DataSet fInitialBField;
	H5::DataSet fSpectralIndex;
	H5::DataSet fL400;
	H5::DataSet fS400;
	H5::DataSet fS1400;
	H5::DataSet fDM;
	H5::DataSet fPositionX;
	H5::DataSet fPositionY;
	H5::DataSet fPositionZ;
	H5::DataSet fVelocityX;
	H5::DataSet fVelocityY;
	H5::DataSet fVelocityZ;
	H5::DataSet fInitialPositionX;
	H5::DataSet fInitialPositionY;
	H5::DataSet fInitialPositionZ;
	H5::DataSet fInitialKickX;
	H5::DataSet fInitialKickY;
	H5::DataSet fInitialKickZ;
	H5::DataSet fGalacticPosB;
	H5::DataSet fGalacticPosL;
	H5::DataSet fGalacticPosDist;
	H5::DataSet fGalacticPosDegB;
	H5::DataSet fGalacticPosDegL;
	H5::DataSet fBeamingFraction;
	H5::DataSet bVisibleInParkes;
	H5::DataSet bVisibleInSKA;
	H5::DataSet bVisibleInSKALC;
};


#if 0
#define BASIC_STAR_SIZE 176
union UBasicStar
{
	SBasicStar Data;
	char  Word[BASIC_STAR_SIZE];

};
#endif
	

class CBasicStar: public SBasicStar
{
public:

	//CBasicStar(){};//TODO

	void PrintStarInfo(std::ostream &out) const;
	void PrintDescriptiveStarInfo(std::ostream &out) const;
	void ResetToInit(const bool bPhysics, const bool bDynamics);///<incomplete should be used with great care. Some variables are random.
	
};

