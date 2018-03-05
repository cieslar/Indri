#pragma once
#include<Definitions.hpp>


class CDMContainer
{
protected:

	//CConfiguration *m_pCfgPtr;
	long int m_nSide;///<Resoultion in the HEALPix scheme
	long int m_nDistNum;///<Number of distance slices
	long int m_nPix;///<number of pixel per distance slice
	float **m_ppDM;///<continous data holder for the Dispersion Measure
	RealType m_fMaxDist;
	RealType m_fDDist;

	void m_PrepareContainer(const long int nSide, const long int nDistNum, const RealType fMaxDist);
public:
	CDMContainer(const long int nSide, const long int nDistNum, const RealType fMaxDist);
	CDMContainer(const char* FileName);
	~CDMContainer();

	RealType GetDM(const RealType fLDeg, const RealType fBDeg, const RealType fDist) const;

	void ReadBinaryGeometry(const char* FileName);
	void ReadBinary(const char* FileName);
	void WriteBinary(const char* FileName);


	void ComputeTheData();
	void ComputeSlice(const long int nSlice);
	void WriteSliceBinary(const long int nSlice, const char* FileName);
	void ReadSliceBinary(const long int nSlice, const char* FileName);
	//RealType FindMaxDiff();
};



RealType GalLDeg2HEALPixPhi(const RealType fLDeg);
RealType HEALPixPhi2GalLDeg(const RealType fPhi);
RealType GalBDeg2HEALPixTheta(const RealType fBDeg);
RealType HEALPixTheta2GalBDeg(const RealType fTheta);

