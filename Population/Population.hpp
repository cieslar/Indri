#pragma once

#include<ConfigContainer.hpp>
#include<PulsarDynamics.hpp>
#include<PulsarPhysics.hpp>

#include<Projections.hpp>
#include<Definitions.hpp>

#include<NE2001.hpp>
#include<Pulsar.hpp>
#include<RadioDetector.hpp>

class CHDFPopulationBuffer
{
protected:
	CConfiguration *m_pCfgPtr;

	static const hsize_t cnDims = 1;
	hsize_t *pOffset;
	hsize_t *pDims;///<buffer size
	hsize_t *pDimsMax;///<max size of dataset

	short *pbBuf;///<"bool" buffer
	double *pfBuf;///<RealType buffer

	H5::H5File *pFile;
	SPopulationDataSet Set;
	CRadioDetector m_RadioObs;
	void CreateDatasets();
public:
	CHDFPopulationBuffer(const char *FileName);
	~CHDFPopulationBuffer();

	void WriteChunk(const CPulsar *pPulsar, const int nIteration);

};


class CPopulation
{
protected:
	RealType m_fMaxTimeStep;

	CConfiguration   *m_pCfgPtr;
	CRandomizer      *m_pRanPtr;
	//CProjection      *m_pProj;

	int      m_nNum;
	RealType m_fMaxAge;
	RealType m_fSunAngle;

	RealType m_PulsarStartingAgeDistribution() const;
    RealType m_PulsarStartingAgeDistribution(const RealType fMinAge, const RealType fMaxAge) const;

	CRadioDetector m_RadioObs;

	bool bHDFIO;
	CHDFPopulationBuffer *pHDFBuf;
public:
	CPulsar	*pPulsar;
		


	CPopulation();
	~CPopulation();
	void CheckVisibleRadio();

	void Init1Per100(const bool bSmearAge=false,const RealType fSmearYearsSigma=5.);
	void Init();
    void Init(const RealType fMinAge, const RealType fMaxAge);

	void Evolve();	
	//void Evolve(const RealType fMaxTimeStep);	
	void Project(const char *FileName);
	void Position(const char *FileName);
	void Write(const char *Filename) const;
	void WriteBinary(const char *Filename) const;
	void ReadBinary(const char *Filename);
//	void ReadBinary(const char *Filename,CConfiguration* pCfgPtr);
	void Zero();
	void ReadBinaryRelevantGeometry(const char *Filename); 
	void ReadHDF5(const char *FileName);
	void WriteHDF5(const char *FileName, const int nIteration);
	void WriteEvolutionTracksHDF5(const char *FileName, const int nStepsPerTrack,const int *pNumery, const int nNumNum);
	void CreateHDF5File(const char *FileName) const;
	void WritePartialHDF5(const char *FileName) const;

	RealType DistanceToObservations();

	void WriteTxtDM(const char *FileName) const;
	void ReadTxtDM(const char *FileName);

	void WriteLost(const char* FileName);
};



class CTestPopulation: public CPopulation
{
public:
	void EnergyTest();
	void WriteReadTest();

	CTestPopulation(){};
	~CTestPopulation(){};

};




