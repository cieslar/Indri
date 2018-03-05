#pragma once

#include<Galaxy.hpp>
#include<ConfigContainer.hpp>
#include<Definitions.hpp>
#include<BasicStar.hpp>
#include<Coordinates.hpp>

class CBinary
{
public:
	int nId;
	
	RealType fTimeBoring;
	RealType fTime1;
	RealType fTime2;
	
	Cartesian Kick1;
	Cartesian Kick2;
 
	Cartesian Pos;
	Cartesian Vel;

	Galactic GalacticPos;
	Galactic GalacticPosDeg;

	CBinary()
	{
		fTimeBoring=fTime1=fTime2=0.;
		nId=0;
	};
	~CBinary(){};
};


class CNSNSBinaryPopulation: protected CGalaxy
{
protected:
	CConfiguration *m_pCfgPtr;
	CRandomizer *m_pRanPtr;
	int m_nNumOfBinaries;
	CBinary *m_pBinary;
	void m_EvolvePositionVerlet(CBinary *pBinary, const RealType fTimeStep) ;
	void m_Kick(CBinary *pBinary, const Cartesian Kick) ;
	void m_EvolveDynamics(CBinary *pBinary, const RealType fEvoTime) ;

public:
	CNSNSBinaryPopulation();
	~CNSNSBinaryPopulation();

	void Read(const char *Filename);
	void Write(const char *Filename);

	void Init();
	void Evolve();
	void Project();

	void WriteProjectionMtx();

};



