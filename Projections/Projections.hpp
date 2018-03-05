#pragma once

#include<iostream>
#include<fstream>
#include<Statistics.hpp>
#include<cmath>
#include<Debug.hpp>
#include<cstring>
#include<Randomizer.hpp>

using namespace std;


//FIXME przeniesc do Tools
#include<limits>
//inline bool DoubleEq(const double fLeft, const double fRight)
bool DoubleEq(const double fLeft, const double fRight);

//FIXME przeniesc do DataStructures
class CDoubleArray2D
{
public:
    size_t nMemMax;
    int nNum;	
    int nWidth;
    int nHight;
    double *pData;
	
    double& Data(const int nX, const int nY);

    CDoubleArray2D(const int nArrayWidth, const int nArrayHight);
    ~CDoubleArray2D();

    void BinWrite(const char *FileName) const;
    void BinRead(const char *FileName);

	void Write(const char *FileName) const;
    void Read(const char *FileName);

    void Zero();
    size_t ComputeMem();
};


class CProjection 
{
protected:
	int m_nResX;
	int m_nResY;
	double m_fZeroX;
	double m_fZeroY;
	double m_fDeltaX;
	double m_fDeltaY;

	CDoubleArray2D m_SquareProjection;

public:
    CProjection(const int nResX, const double fLZero=0.);
    ~CProjection(){};

    void Write(const char *FileName);
    void BinWrite(const char *FileName);
    void Zero();
    void AddPoint(const double fL, const double fB);
	void SphereDensify();
	void Norm();

};

class CMollweideProjection : public CProjection
{
protected:
	CDoubleArray2D m_MollweideProjection;
	double m_fRadius;
	double m_ComputeTheta(const double fB);
double m_ComputeThetaBisec(const double fB);
public:
    	CMollweideProjection(const int nResX, const double fLZero=0.);
    	~CMollweideProjection(){};

    	void Write(const char *FileName);
    	void BinWrite(const char *FileName);
    	void Zero();
	void Project();
};


