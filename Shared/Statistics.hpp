#pragma once
#include<Definitions.hpp>
#include<ModelParameters.hpp>
#include<CatalogueMasks.hpp>

//FIXME jesli bede robil refaktoryzacje to CSpaceDensity i CRadioDensity można zwinąć do jednego z kontruktorem od struktury opisującej parametry wykresu.
//TODO za każdym razem jak zmienię jedną z nich, muszą pamiętać aby zmienić też drugą


class CDensityVector
{
protected:
	int nTotalNum;
	CConfiguration *m_pCfgPtr;
	void AddPeriod(const RealType Period);
	void AddPeriodDot(const RealType PeriodDot);
	void AddS1400(const RealType S1400);
	RealType fDLogP;
	RealType fDLogPdot;
	RealType fDS1400;

	RealType fLogPLeft;
	RealType fLogPRight;
	RealType fLogPdotLeft;
	RealType fLogPdotRight;
	RealType fLogS1400Left;
	RealType fLogS1400Right;
public:

	RealType fPeriodNumFraction;
	RealType fPeriodDotNumFraction;
	RealType fS1400NumFraction;


	int nPeriodSize;
	int nPeriodDotSize;
	int nS1400Size;
	RealType *pPeriodLog;
	RealType *pPeriodDotLog;
	RealType *pS1400Log;
	
	RealType *pPeriodLogCD;
	RealType *pPeriodDotLogCD;
	RealType *pS1400LogCD;

	void Norm(const int nNum);
	void Norm();
	void Zero();
	void Cumulate();

	void AddPoint(const CPulsar PSR);
	void AddPoint(const SCatalogueEntry OBS);
	CDensityVector();
	~CDensityVector();
};

class CSpaceDensity
{
protected:
	int m_nResGalL;
	int m_nResGalKsi;///sin(M_PI+GalKsi)
	int m_nResDM;

	RealType m_fZeroGalL;
	RealType m_fZeroGalKsi;
	RealType m_fZeroDM;

	RealType m_fMaxGalL;
	RealType m_fMaxGalKsi;
	RealType m_fMaxDM;

	RealType m_fDGalL;
	RealType m_fDGalKsi;
	RealType m_fDDM;

	RealType *m_pfDensity;

public:
	int m_nAddedPSRs;
	int m_nTotalPSRs;
	RealType fDensityElement;
	void AddPoint(const CPulsar PSR);
	void AddPoint(const RealType fGalL, const RealType fGalKsi, const RealType fDM, const RealType fProbability);
	RealType GetDensity(const RealType fGalL, const RealType fGalKsi, const RealType fDM);
	RealType GetDensity(const int nGalL, const int nGalKsi, const int nDM);
	~CSpaceDensity();
	CSpaceDensity(const RealType fMaxDM, const unsigned int nResGalL, const unsigned int nResGalKsi, const unsigned int nResDM );

	void Norm();
	void Densify();
	void Zero();

	void Write(const char *FileName);

	RealType Check();

};





class CRadioDensity
{
protected:
	int m_nResLogP;
	int m_nResLogPdot;
	int m_nResLogS1400;

	RealType m_fZeroLogP;
	RealType m_fZeroLogPdot;
	RealType m_fZeroLogS1400;

	RealType m_fMaxLogP;
	RealType m_fMaxLogPdot;
	RealType m_fMaxLogS1400;

	RealType m_fDLogP;
	RealType m_fDLogPdot;
	RealType m_fDLogS1400;

	RealType *m_pfDensity;


public:
	int m_nAddedPSRs;
	int m_nTotalPSRs;
	RealType fDensityElement;

	void AddValueAtPoint(const int nLogP, const int nLogPdot, const int nLogS1400, const RealType fValue);

	void AddPoint(const CPulsar PSR);
	void AddPoint(const RealType fLogP, const RealType fLogPdot, const RealType fLogS1400, const RealType fProbability);
	RealType GetDensity(const RealType fLogP, const RealType fLogPdot, const RealType fLogS1400);
	RealType GetDensity(const int nLogP, const int nLogPdot, const int nLogS1400);
	~CRadioDensity();
	CRadioDensity(	const RealType fZeroLogP, const RealType fZeroLogPdot, const RealType fZeroLogS1400, const RealType fMaxLogP, const RealType fMaxLogPdot, const RealType fMaxLogS1400, const unsigned int nResLogP, const unsigned int nResLogPdot, const unsigned int nResLogS1400 );

	void Norm();
	void Densify();
	void Zero();
	void Write(const char *FileName);
	RealType Check();

};




class CDensity
{
public:
	virtual void Norm();//1/N
	virtual void Densify();//physical dV element 
	virtual void Zero();
	virtual RealType GetDensitySpace(const SProbabilitySpacePointSpace SpacePoint);
	virtual RealType GetDensityRadio(const SProbabilitySpacePointRadio SpacePoint);

	virtual RealType GetDensity(const SProbabilitySpacePointRadio SpacePointRadio, const SProbabilitySpacePointSpace SpacePointSpace); 
	virtual void AddPoint(const CPulsar Pulsar);
	virtual void AddValueAtPoint(const SProbabilitySpacePointRadio SpacePointRadio, const RealType fValue);

	virtual void Write(const char* FilePrefix);
	virtual void Check();//check the normalization

	virtual RealType GetEfficiency();

	
	CDensity();
	~CDensity();
};


class CDensity3x1D : public CDensity
{
protected:
	CConfiguration *m_pCfgPtr;
	//CDistribution *m_ppDistr;
public:
	CDensityVector Vector;
	void Norm();
	void Densify();
	void Zero();
	void AddPoint(const CPulsar Pulsar);
	CDensity3x1D();
	~CDensity3x1D();
};

class CDensity2x3D : public CDensity
{
protected:
	CConfiguration *m_pCfgPtr;
	CRadioDensity *m_pRadio;
	CSpaceDensity *m_pSpace;
public:
	RealType GetEfficiency();
	void Norm();
	void Densify();
	void Zero();
	RealType GetDensitySpace(const SProbabilitySpacePointSpace SpacePoint);
	RealType GetDensityRadio(const SProbabilitySpacePointRadio SpacePoint);

	RealType GetDensity(const SProbabilitySpacePointRadio SpacePointRadio, const SProbabilitySpacePointSpace SpacePointSpace); 
	void AddPoint(const CPulsar Pulsar);
	void AddValueAtPoint(const SProbabilitySpacePointRadio SpacePointRadio, const RealType fValue);
	void Write(const char* FilePrefix);
	void Check();
	
	CDensity2x3D();
	~CDensity2x3D();
};

class CDensity2x3DSparse : public CDensity
{
protected:

	int m_nAddedPSRsRadio;
	int m_nTotalPSRsRadio;
	int m_nAddedPSRsSpace;
	int m_nTotalPSRsSpace;
	
	RealType m_fVolumeElementRadio;
	RealType m_fVolumeElementSpace;

	RealType m_ComputeVolumeElementRadio() const;
	RealType m_ComputeVolumeElementSpace() const;
	void m_AddToMap(const SProbabilitySpacePointRadio Point);
	void m_AddToMap(const SProbabilitySpacePointSpace Point);


	CConfiguration *m_pCfgPtr;

public:
	RealType GetEfficiency();
	RealType fDensityElement;

	std::map<SProbabilitySpacePointRadio, RealType> MapRadio;
	std::map<SProbabilitySpacePointSpace, RealType> MapSpace;

	void Zero();

	void CopyMapKeys(const CObservations* pObservations);
	RealType GetDensitySpace(const SProbabilitySpacePointSpace SpacePoint);
	RealType GetDensityRadio(const SProbabilitySpacePointRadio SpacePoint);

	RealType GetDensity(const SProbabilitySpacePointRadio SpacePointRadio, const SProbabilitySpacePointSpace SpacePointSpace); 

	void AddPoint(const CPulsar Pulsar);
	
	CDensity2x3DSparse(const CObservations* pObservations);
	~CDensity2x3DSparse();

	//Empty functions
	void Check();
	void Norm();
	void Densify();
	void Write(const char* FilePrefix);
};

class CDensity4DSparse : public CDensity
{
protected:

	int m_nAddedPSRs;
	int m_nTotalPSRs;
	
	RealType m_fVolumeElement;

	RealType m_ComputeVolumeElement() const;
	void m_AddToMap(const SProbabilitySpacePoint4D Point);


	CConfiguration *m_pCfgPtr;

public:
	RealType GetEfficiency();
	RealType fDensityElement;

	std::map<SProbabilitySpacePoint4D, RealType> Map4D;

	void Zero();

	void CopyMapKeys(const CObservations* pObservations);
	RealType GetDensity(const SProbabilitySpacePoint4D SpacePoint);

	void AddPoint(const CPulsar Pulsar);
	
	CDensity4DSparse(const CObservations* pObservations);
	~CDensity4DSparse();

	//Empty functions
	void Check();
	void Norm();
	void Densify();
	void Write(const char* FilePrefix);
};



struct SCumulDensity
{
	int nX;
	int nY;
	double fValue; 
};

bool operator<(const SCumulDensity &A, const SCumulDensity &B);
//FIXME do refaktoryzacji (CS)
class CAreaDensity
{
public:
	int nResX;
	int nResY;
	double fZeroX;
	double fZeroY;
	double fMaxX;
	double fMaxY;
	double fDX;
	double fDY;
	
	double **fDensity;//Musze pamietac, ze jestem ograniczony 16 cyframi znaczacymi przy dodawaniu +1 dla doubli
	//fDensity[y][x];

	unsigned char *pSigma;
	SCumulDensity *pCumulative;
	bool bCumulative;

	void Sigmify();
	void WriteSigma(const char* Filename);
	unsigned char Sigma(const double fX, const double fY);

	CAreaDensity(const double ZeroX, const double ZeroY, const double MaxX, const double MaxY, const unsigned int ResX, const unsigned int ResY);
	~CAreaDensity();
	void AddPoint(const double fX, const double fY);

	void Densify();//Dzieli przez element powierzchni
	void Norm();//Normuje do 1
	void Write(const char *FileName);
	//WriteLog();

	void Reset();
};

class CCubicDensity
{
public:
	int nResX;
	int nResY;
	int nResZ;
	double fZeroX;
	double fZeroY;
	double fZeroZ;
	double fMaxX;
	double fMaxY;
	double fMaxZ;
	double fDX;
	double fDY;
	double fDZ;
	
	double *pfDensity;

	

	CCubicDensity(const double ZeroX, const double ZeroY, const double ZeroZ, const double MaxX, const double MaxY, const double MaxZ, const  int ResX, const  int ResY, const  int ResZ);
	~CCubicDensity();
	void AddPoint(const double fX, const double fY, const double fZ);

	int Index(const double fX, const double fY, const double fZ);
	void Zero();
	void Densify();//Dzieli przez element objetosci
	void Norm();//Normuje do 1
	void Write(const char *FileName);
};




class CHistogram
{
private:
	int 	*m_pCount;
	double 	*m_pDensity;
	double 	m_fMinValue;
	double 	m_fMaxValue;
	double 	m_fBinWidth;
	int	m_nBinNum;

public:
	CHistogram(const double fMinValue, const double fMaxValue, const double fBinWidth);
	CHistogram(const double fMinValue, const double fMaxValue, const int nBins);

	~CHistogram();

	void AddPoint(const double fValue);

	void Zero();
	
	void Densify();//Dzieli przez szerokosc binu
	void Norm();
	void Norm(const double fNormValue);
	void OnlyCount();

	void Write(const char *FileName);
	
};



class CCumulativeDistribution
{
private:
	int 	*m_pCount;
	double 	*m_pCumuDist;

	double 	m_fMinValue;
	double 	m_fMaxValue;
	double 	m_fBinWidth;
	int	m_nBinNum;

public:
	CCumulativeDistribution(const double fMinValue, const double fMaxValue, const int nResolution);
	~CCumulativeDistribution();

	void AddPoint(const double fValue);

	void Zero();
	
	void Norm();
	void Norm(const double fNormValue);
	void Cumulify();

	void Write(const char *FileName);
	
};



RealType DValue3D(const CDensityVector *pVA, const CDensityVector *pVB);


	

