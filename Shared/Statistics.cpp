#include<Statistics.hpp>

#include<fstream>
#include<ConfigContainer.hpp>
#include<Pulsar.hpp>


CHistogram::CHistogram(const double fMinValue, const double fMaxValue, const double fBinWidth)
{
	m_fMinValue = fMinValue;
	m_fMaxValue = fMaxValue;
	m_fBinWidth = fBinWidth;
	m_nBinNum   = static_cast<int>( (m_fMaxValue-m_fMinValue) / fBinWidth ) + 1;

	m_pCount   = new int[m_nBinNum];
	m_pDensity = new double[m_nBinNum];

	Zero();
}

CHistogram::CHistogram(const double fMinValue, const double fMaxValue, const int nBins)
{
	m_fMinValue = fMinValue;
	m_fMaxValue = fMaxValue;
	m_nBinNum   = nBins;
	m_fBinWidth = (m_fMaxValue-m_fMinValue) / (m_nBinNum - 1);

	m_pCount   = new int[m_nBinNum];
	m_pDensity = new double[m_nBinNum];

	Zero();
}

CHistogram::~CHistogram()
{
	delete [] m_pCount;
	delete [] m_pDensity;
}

void CHistogram::AddPoint(const double fValue)
{
	if( (fValue <= m_fMaxValue) && (fValue >= m_fMinValue) )
	{
		int index = static_cast<int>( (fValue-m_fMinValue) / m_fBinWidth );
		m_pCount[index]++;
	}
}

void CHistogram::Zero()
{
	for(int i(0);i<m_nBinNum;i++)
	{
		m_pCount[i]=0;
		m_pDensity[i]=0.;		
	}
}

void CHistogram::OnlyCount()
{
	for(int i(0);i<m_nBinNum;i++)
	{
		m_pDensity[i]=m_pCount[i];
	}
}

void CHistogram::Densify()
{
	for(int i(0);i<m_nBinNum;i++)
	{
		m_pDensity[i]=m_pCount[i]/m_fBinWidth;
	}
}

void CHistogram::Norm()
{
	double fTotalValue(0.);
	for(int i(0);i<m_nBinNum;i++) fTotalValue+=m_pDensity[i];
	Norm(fTotalValue);
}

void CHistogram::Norm(const double fNormValue)
{
	for(int i(0);i<m_nBinNum;i++)
	{
		m_pDensity[i]=m_pCount[i]/fNormValue;
	}
}

void CHistogram::Write(const char *FileName)
{
	std::ofstream out(FileName);
	for(int i(0);i<m_nBinNum;i++)
	{
		out<<m_fMinValue + (0.5 + i)*m_fBinWidth <<" "<<m_pDensity[i]<<std::endl;
	}

	out.close();

}

//////////////////////////////////////////////////
//////////////////////////////////////////////////



CCumulativeDistribution::CCumulativeDistribution(const double fMinValue, const double fMaxValue, const int nResolution) 
{
	m_fMinValue = fMinValue;
	m_fMaxValue = fMaxValue;
	m_nBinNum   = nResolution;
	m_fBinWidth = (m_fMaxValue-m_fMinValue) / (m_nBinNum-1);//bo jest zero i ostatni punkt, rozdzielczosc tyczy sie ilosci przedzialow.

	m_pCount    = new int[m_nBinNum];
	m_pCumuDist = new double[m_nBinNum];

	Zero();
}

CCumulativeDistribution::~CCumulativeDistribution()
{
	delete [] m_pCount;
	delete [] m_pCumuDist;
}

void CCumulativeDistribution::AddPoint(const double fValue)
{
	if( (fValue <= m_fMaxValue) && (fValue >= m_fMinValue) )
	{
		int index = static_cast<int>( (fValue-m_fMinValue) / m_fBinWidth );
		m_pCount[index]++;
	}
}

void CCumulativeDistribution::Zero()
{
	for(int i(0);i<m_nBinNum;i++)
	{
		m_pCount[i]=0;
		m_pCumuDist[i]=0.;		
	}
}

void CCumulativeDistribution::Cumulify()
{
	m_pCumuDist[0]=m_pCount[0];
	for(int i(1);i<m_nBinNum;i++)
	{	
		m_pCumuDist[i]=m_pCount[i]+m_pCumuDist[i-1];
	}
}

void CCumulativeDistribution::Norm()
{
	Norm(m_pCumuDist[m_nBinNum-1]);
}

void CCumulativeDistribution::Norm(const double fNormValue)
{
	for(int i(0);i<m_nBinNum;i++)
	{
		m_pCumuDist[i]=m_pCumuDist[i]/fNormValue;
	}
}

void CCumulativeDistribution::Write(const char *FileName)
{
	std::ofstream out(FileName);
	for(int i(0);i<m_nBinNum;i++)
	{
		out<<m_fMinValue + i*m_fBinWidth <<" "<<m_pCumuDist[i]<<std::endl;
	}

	out.close();

}




////////////////////////////////////////////////
////////////////////////////////////////////////
#include<algorithm>
#include <iterator>
bool operator<(const SCumulDensity &A, const SCumulDensity &B)
{
	if(A.fValue<B.fValue) return true;
	return false;
}

void CAreaDensity::Sigmify()
{
	bCumulative = true;
	pCumulative = new SCumulDensity[nResX*nResY];
	pSigma = new unsigned char[nResX*nResY];
//znajdz maksimum

	for(int X(0);X<nResX;X++) for(int Y(0);Y<nResY;Y++)
	{
		pCumulative[X+nResX*Y].nX=X;
		pCumulative[X+nResX*Y].nY=Y;
		pCumulative[X+nResX*Y].fValue=fDensity[Y][X];
	}

	sort(pCumulative, pCumulative + nResX*nResY );

	const double f3Sigma(0.9973);
	const double f2Sigma(0.9545);
	const double f1Sigma(0.6827);
	
//TODO sprawdzic jakie sa wartosci dla przedzialow ufnosci!

	int nNum=nResX*nResY;	
	for(int i(nNum-1);i>=0;i--)
	{
		if(i<nNum-1) pCumulative[i].fValue+=pCumulative[i+1].fValue;

		pSigma[pCumulative[i].nX + nResX*pCumulative[i].nY]=0;
		if(pCumulative[i].fValue<=f3Sigma)
		{
			pSigma[pCumulative[i].nX + nResX*pCumulative[i].nY]=3;
			if(pCumulative[i].fValue<=f2Sigma)
			{
				pSigma[pCumulative[i].nX + nResX*pCumulative[i].nY]=2;
				if(pCumulative[i].fValue<=f1Sigma)
				{
					pSigma[pCumulative[i].nX + nResX*pCumulative[i].nY]=1;
				}
			}
			
		}
	}

}

void CAreaDensity::WriteSigma(const char* Filename)
{	
	std::ofstream out(Filename);
	for(unsigned int j(0);j<nResY;j++)
	{
		for(unsigned int i(0);i<nResX;i++) 
		{
			out<<static_cast<unsigned int>(pSigma[i+j*nResX])<<" ";
		}
		out<<std::endl;
	}
	out.close();
}



unsigned char CAreaDensity::Sigma(const double fX, const double fY)
{
	int X=static_cast<int>((fX-fZeroX)/fDX);
	int Y=static_cast<int>((fY-fZeroY)/fDY);
	if(X>=0 && Y>=0 && Y<nResY && X<nResX && bCumulative) return pSigma[X+nResX*Y];
	return 0;
}



CAreaDensity::CAreaDensity(const double ZeroX, const double ZeroY, const double MaxX, const double MaxY, const unsigned int ResX, const unsigned int ResY)
{
	fZeroX=ZeroX;
	fMaxX=MaxX;
	nResX=ResX;
	fZeroY=ZeroY;
	fMaxY=MaxY;
	nResY=ResY;
	fDX=(fMaxX-fZeroX)/nResX;
	fDY=(fMaxY-fZeroY)/nResY;
	bCumulative=false;

	fDensity=new double*[nResY];//FIXME raz, ze codestyle a dwa na accessor: double& Data(x,y);
	for(unsigned int i(0);i<nResY;i++) fDensity[i]=new double[nResX];
	for(unsigned int i(0);i<nResY;i++) for(unsigned int j(0);j<nResX;j++) fDensity[i][j]=0.;
	
}


CAreaDensity::~CAreaDensity()
{
	for(unsigned int i(0);i<nResY;i++) delete [] fDensity[i];
	delete [] fDensity;
	if(bCumulative)
	{
		delete [] pCumulative;
		delete [] pSigma;
	}
}

void CAreaDensity::AddPoint(const double fX, const double fY)
{
	//FIXME dodac checkera!
	int X=static_cast<int>((fX-fZeroX)/fDX);
	int Y=static_cast<int>((fY-fZeroY)/fDY);
	if(X>=0 && Y>=0 && Y<nResY && X<nResX) fDensity[Y][X]++;
}

void CAreaDensity::Densify()
{
	for(unsigned int i(0);i<nResY;i++) for(unsigned int j(0);j<nResX;j++) fDensity[i][j]/=fDX*fDY;
}

void CAreaDensity::Norm()
{
	double fAll(0.);//technicznie bede robil blad dla duzych zbiorow ze wzgledu na kolejnosc dodawania (przekroczenie zakresu double)
	for(unsigned int i(0);i<nResY;i++) for(unsigned int j(0);j<nResX;j++) fAll+=fDensity[i][j];
	for(unsigned int i(0);i<nResY;i++) for(unsigned int j(0);j<nResX;j++) fDensity[i][j]/=fAll;
}
void CAreaDensity::Reset()
{
	for(unsigned int i(0);i<nResY;i++) 
	{
		for(unsigned int j(0);j<nResX;j++)
		{
			fDensity[i][j]=0.;
		}
	}
}
void CAreaDensity::Write(const char *FileName)
{
	std::ofstream out(FileName);
	for(unsigned int i(0);i<nResY;i++) 
	{
		for(unsigned int j(0);j<nResX;j++)
		{
			out<<fDensity[i][j]<<" ";
		}
		out<<std::endl;
	}
	out.close();
}


///////////////////////////////////////////////////////////////////
int CCubicDensity::Index(const double fX, const double fY, const double fZ)
{
	int X=static_cast<int>((fX-fZeroX)/fDX);
	int Y=static_cast<int>((fY-fZeroY)/fDY);
	int Z=static_cast<int>((fZ-fZeroZ)/fDZ);
	int index = X + Y*nResX + Z*nResX*nResY ;
	if(index<nResX*nResY*nResZ) return index;
	return -1;
}

void CCubicDensity::Zero()
{
	for(int i(0);i<nResX*nResY*nResZ;i++) pfDensity[i]=0.;

}


CCubicDensity::CCubicDensity(const double ZeroX, const double ZeroY, const double ZeroZ, const double MaxX, const double MaxY, const double MaxZ, const  int ResX, const  int ResY, const  int ResZ)
{
	fZeroX=ZeroX;
	fMaxX=MaxX;
	nResX=ResX;
	fZeroY=ZeroY;
	fMaxY=MaxY;
	nResY=ResY;
	fZeroZ=ZeroZ;
	fMaxZ=MaxZ;
	nResZ=ResZ;
	fDZ=(fMaxZ-fZeroZ)/nResZ;
	fDX=(fMaxX-fZeroX)/nResX;
	fDY=(fMaxY-fZeroY)/nResY;

	pfDensity=new double[nResY*nResX*nResZ];

	Zero();
}


CCubicDensity::~CCubicDensity()
{
	delete [] pfDensity;
}

void CCubicDensity::AddPoint(const double fX, const double fY, const double fZ)
{
	int i=Index(fX,fY,fZ);
	if(i>=0) pfDensity[i]++;
}

void CCubicDensity::Densify()
{
	for(int i(0);i<nResX;i++) for(int j(0);j<nResY;j++) for(int k(0);k<nResZ;k++) pfDensity[i +j*nResX +k*nResX*nResY]/=fDX*fDY*fDZ;
}

void CCubicDensity::Norm()
{
	double fAll(0.);//technicznie bede robil blad dla duzych zbiorow ze wzgledu na kolejnosc dodawania (przekroczenie zakresu double)
	for(int i(0);i<nResX;i++) for(int j(0);j<nResY;j++) for(int k(0);k<nResZ;k++) fAll+=pfDensity[i +j*nResX +k*nResX*nResY];
	for(int i(0);i<nResX;i++) for(int j(0);j<nResY;j++) for(int k(0);k<nResZ;k++) pfDensity[i +j*nResX +k*nResX*nResY]/=fAll;
}

void CCubicDensity::Write(const char *FileName)
{
	std::ofstream out(FileName);
	for(int i(0);i<nResX;i++) 
	{
		for(int j(0);j<nResY;j++)
		{
			for(int k(0);k<nResZ;k++)
			{
				if(pfDensity[i +j*nResX +k*nResX*nResY]>0.) out<<i*fDX<<" "<<j*fDY<<" "<<k*fDZ<<" "<<pfDensity[i +j*nResX +k*nResX*nResY]<<endl;
			}
		}
	}
	out.close();
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//CRadioDensity


CRadioDensity::~CRadioDensity()
{
	delete [] m_pfDensity;
}

CRadioDensity::CRadioDensity(	const RealType fZeroLogP, const RealType fZeroLogPdot, const RealType fZeroLogS1400, const RealType fMaxLogP, const RealType fMaxLogPdot, const RealType fMaxLogS1400, const unsigned int nResLogP, const unsigned int nResLogPdot, const unsigned int nResLogS1400 )
{
	m_nResLogP = nResLogP;
	m_nResLogPdot = nResLogPdot;
	m_nResLogS1400 = nResLogS1400;

	m_fZeroLogP = fZeroLogP;
	m_fZeroLogPdot = fZeroLogPdot;
	m_fZeroLogS1400 = fZeroLogS1400;

	m_fMaxLogP = fMaxLogP;
	m_fMaxLogPdot = fMaxLogPdot;
	m_fMaxLogS1400 = fMaxLogS1400;

	m_fDLogP = (fMaxLogP - fZeroLogP)/nResLogP;
	m_fDLogPdot = (fMaxLogPdot - fZeroLogPdot)/nResLogPdot;
	m_fDLogS1400 = (fMaxLogS1400 - fZeroLogS1400)/nResLogS1400;

	m_pfDensity = new RealType[nResLogP * nResLogPdot * nResLogS1400];

	Zero();

}

void CRadioDensity::Zero()
{
	fDensityElement=1.;
	m_nTotalPSRs=0;
	m_nAddedPSRs=0;
	for(int i(0);i<m_nResLogP * m_nResLogPdot * m_nResLogS1400;i++) m_pfDensity[i]=0.;
}


void CRadioDensity::Norm()
{
    double fRes(0.);
	for(int i(0);i<m_nResLogP * m_nResLogPdot * m_nResLogS1400;i++) fRes+=m_pfDensity[i];
	if(fRes>0.) for(int i(0);i<m_nResLogP * m_nResLogPdot * m_nResLogS1400;i++) m_pfDensity[i]/=fRes;
	//if(m_nTotalPSRs > 0) fDensityElement /= m_nTotalPSRs;
}

void CRadioDensity::Densify()
{
	RealType fDQubic = m_fDLogP * m_fDLogPdot * m_fDLogS1400;
	fDensityElement /= fDQubic;
}

RealType CRadioDensity::Check()
{
	RealType fTotal(0.);
	for(int i(0);i<m_nResLogP * m_nResLogPdot * m_nResLogS1400;i++) fTotal+=m_pfDensity[i];
	return fTotal;

}
void CRadioDensity::AddPoint(const RealType fLogP, const RealType fLogPdot, const RealType fLogS1400, const RealType fProbability)
{

	int X=static_cast<int>((fLogP-m_fZeroLogP)/m_fDLogP);
	int Y=static_cast<int>((fLogPdot-m_fZeroLogPdot)/m_fDLogPdot);
	int Z=static_cast<int>((fLogS1400-m_fZeroLogS1400)/m_fDLogS1400);
	//cout<<"Dodawanie punktu (Radio) f: "<<fLogP<<" "<<fLogPdot<<" "<<fLogS1400<<endl;	
	//cout<<"Dodawanie punktu (Radio) n: "<<X<<" "<<Y<<" "<<Z<<endl;	

	if(X>=0 && Y>=0 && Z>=0 && X<m_nResLogP && Y<m_nResLogPdot && Z<m_nResLogS1400) 
	{
		m_nAddedPSRs++;
		m_pfDensity[X + Y*m_nResLogP + Z*m_nResLogP*m_nResLogPdot]+=fProbability;
		//cout<<" Dodawanie punktu (Radio) f: "<<fLogP<<" "<<fLogPdot<<" "<<fLogS1400<<endl;	
		//cout<<" Dodawanie punktu (Radio) n: "<<X<<" "<<Y<<" "<<Z<<endl;	
		//cout<<" fProbability: "<<fProbability<<endl;
	}
	m_nTotalPSRs++;

}

void CRadioDensity::AddPoint(const CPulsar PSR)
{
	AddPoint(log10(PSR.fPeriod),log10(PSR.fPeriodDot), log10(PSR.fS1400),PSR.fProbability);
}


RealType CRadioDensity::GetDensity(const RealType fLogP, const RealType fLogPdot, const RealType fLogS1400)
{
	int X=static_cast<int>((fLogP-m_fZeroLogP)/m_fDLogP);
	int Y=static_cast<int>((fLogPdot-m_fZeroLogPdot)/m_fDLogPdot);
	int Z=static_cast<int>((fLogS1400-m_fZeroLogS1400)/m_fDLogS1400);
	
	if(X>=0 && Y>=0 && Z>=0 && X<m_nResLogP && Y<m_nResLogPdot && Z<m_nResLogS1400) 
	{
		return m_pfDensity[X + Y*m_nResLogP + Z*m_nResLogP*m_nResLogPdot]*fDensityElement;
	}
	return 0.;

}

void CRadioDensity::AddValueAtPoint(const int nLogP, const int nLogPdot, const int nLogS1400, const RealType fValue)
{
  	if(nLogP>=0 && nLogPdot>=0 && nLogS1400>=0 && nLogP<m_nResLogP && nLogPdot<m_nResLogPdot && nLogS1400<m_nResLogS1400) 
    {
    	m_nAddedPSRs++;
        m_pfDensity[nLogP + nLogPdot*m_nResLogP + nLogS1400*m_nResLogP*m_nResLogPdot]+=fValue;
    }
#if 0
    //This will hapen for SKA.
    else
    {
        cout<<"This should not happen!"<<endl;
        cout<<nLogP<<" "<<nLogPdot*m_nResLogP<<" "<<nLogS1400<<" "<<m_nResLogP*m_nResLogPdot<<endl;
    }
#endif
    m_nTotalPSRs++;
}



RealType CRadioDensity::GetDensity(const int nLogP, const int nLogPdot, const int nLogS1400)
{
	if(nLogP>=0 && nLogPdot>=0 && nLogS1400>=0 && nLogP<m_nResLogP && nLogPdot<m_nResLogPdot && nLogS1400<m_nResLogS1400) 
	{
		return m_pfDensity[nLogP + nLogPdot*m_nResLogP + nLogS1400*m_nResLogP*m_nResLogPdot]*fDensityElement;
	}
	return 0.;

}


void CRadioDensity::Write(const char *FileName)
{
	ofstream out(FileName);

	for(int i(0);i<m_nResLogP;i++) for(int j(0);j<m_nResLogPdot;j++) for(int k(0);k<m_nResLogS1400;k++) 
	if(m_pfDensity[i+j*m_nResLogP+k*m_nResLogP*m_nResLogPdot]>0.)
		out<<m_fZeroLogP +i*m_fDLogP<<" "<<m_fZeroLogPdot+j*m_fDLogPdot<<" "<<m_fZeroLogS1400+k*m_fDLogS1400<<" "<<m_pfDensity[i+j*m_nResLogP+k*m_nResLogP*m_nResLogPdot]<<endl;

	out.close();
}

///////////////////////////////////////////////////////////////////////
CSpaceDensity::~CSpaceDensity()
{
	delete [] m_pfDensity;
}

CSpaceDensity::CSpaceDensity(const RealType fMaxDM, const unsigned int nResGalL, const unsigned int nResGalKsi, const unsigned int nResDM )
{
	m_nResGalL = nResGalL;
	m_nResGalKsi = nResGalKsi;
	m_nResDM = nResDM;

	m_fZeroGalL = 0.;
	m_fZeroGalKsi = -1.;
	m_fZeroDM = 0.;

	m_fMaxGalL = 360.;//FIXME poprawic abym nie mial zdublowanego 0
	m_fMaxGalKsi = 1.;
	m_fMaxDM = fMaxDM;

	m_fDGalL = (m_fMaxGalL - m_fZeroGalL)/nResGalL;
	m_fDGalKsi = (m_fMaxGalKsi - m_fZeroGalKsi)/nResGalKsi;
	m_fDDM = (fMaxDM - m_fZeroDM)/nResDM;

	m_pfDensity = new RealType[nResGalL * nResGalKsi * nResDM];

	Zero();

}

void CSpaceDensity::Zero()
{
	fDensityElement=1.;
	m_nTotalPSRs=0;
	m_nAddedPSRs=0;
	for(int i(0);i<m_nResGalL * m_nResGalKsi * m_nResDM;i++) m_pfDensity[i]=0.;
}


void CSpaceDensity::Norm()
{
    double fRes(0.);
	for(int i(0);i<m_nResGalL * m_nResGalKsi * m_nResDM;i++) fRes+=m_pfDensity[i];
	if(fRes>0.) for(int i(0);i<m_nResGalL * m_nResGalKsi * m_nResDM;i++) m_pfDensity[i]/=fRes;
	//if(m_nTotalPSRs > 0) fDensityElement /= m_nTotalPSRs;
//	for(int i(0);i<m_nResGalL * m_nResGalKsi * m_nResDM;i++) fTotal+=m_pfDensity[i];
//	if(fTotal>0.) for(int i(0);i<m_nResGalL * m_nResGalKsi * m_nResDM;i++) m_pfDensity[i]/=fTotal;

}
void CSpaceDensity::Densify()
{
	RealType fDQubic = m_fDGalL * m_fDGalKsi * m_fDDM;
	fDensityElement /= fDQubic;
}

RealType CSpaceDensity::Check()
{
	RealType fTotal(0.);
	for(int i(0);i<m_nResGalL * m_nResGalKsi * m_nResDM;i++) fTotal+=m_pfDensity[i];

	return fTotal;

}


void CSpaceDensity::AddPoint(const RealType fGalL, const RealType fGalKsi, const RealType fDM, const RealType fProbability)
{

	int X=static_cast<int>((fGalL-m_fZeroGalL)/m_fDGalL);
	int Y=static_cast<int>((fGalKsi-m_fZeroGalKsi)/m_fDGalKsi);
	int Z=static_cast<int>((fDM-m_fZeroDM)/m_fDDM);
	//cout<<"Dodawanie punktu (Space): "<<X<<" "<<Y<<" "<<Z<<endl;	
	
	if(X>=0 && Y>=0 && Z>=0 && X<m_nResGalL && Y<m_nResGalKsi && Z<m_nResDM) 
	{
		m_nAddedPSRs++;
		m_pfDensity[X + Y*m_nResGalL + Z*m_nResGalL*m_nResGalKsi]+=fProbability;
	}
	m_nTotalPSRs++;

}
void CSpaceDensity::AddPoint(const CPulsar PSR)
{
	//cout<<PSR.GalacticPosDeg.fL<<" "<<PSR.GalacticPosDeg.fB<<" "<< PSR.GalacticPos.fL<< " "<<PSR.GalacticPos.fB<<" "<<sin(PSR.GalacticPos.fB)<<endl;
	AddPoint(PSR.GalacticPosDeg.fL, sin(PSR.GalacticPos.fB), PSR.fDM, PSR.fProbability);


}

RealType CSpaceDensity::GetDensity(const RealType fGalL, const RealType fGalKsi, const RealType fDM)
{
	int X=static_cast<int>((fGalL-m_fZeroGalL)/m_fDGalL);
	int Y=static_cast<int>((fGalKsi-m_fZeroGalKsi)/m_fDGalKsi);
	int Z=static_cast<int>((fDM-m_fZeroDM)/m_fDDM);
	
	if(X>=0 && Y>=0 && Z>=0 && X<m_nResGalL && Y<m_nResGalKsi && Z<m_nResDM) 
	{
		return m_pfDensity[X + Y*m_nResGalL + Z*m_nResGalL*m_nResGalKsi]*fDensityElement;
	}
	return 0.;

}


RealType CSpaceDensity::GetDensity(const int nGalL, const int nGalKsi, const int nDM)
{
	if(nGalL>=0 && nGalKsi>=0 && nDM>=0 && nGalL<m_nResGalL && nGalKsi<m_nResGalKsi && nDM<m_nResDM) 
	{
		return m_pfDensity[nGalL + nGalKsi*m_nResGalL + nDM*m_nResGalL*m_nResGalKsi]*fDensityElement;
	}
	return 0.;

}


void CSpaceDensity::Write(const char *FileName)
{
	ofstream out(FileName);

	for(int i(0);i<m_nResGalL;i++)	for(int j(0);j<m_nResGalKsi;j++) for(int k(0);k<m_nResDM;k++) 
	if(m_pfDensity[i+j*m_nResGalL+k*m_nResGalL*m_nResGalKsi]>0.)
	out<< ((m_fZeroGalL + i*m_fDGalL)>180.? (m_fZeroGalL + i*m_fDGalL)-360. : (m_fZeroGalL + i*m_fDGalL)) <<" "<<m_fZeroGalKsi + j*m_fDGalKsi<<" "<<m_fZeroDM + k*m_fDDM<<" "<<m_pfDensity[i+j*m_nResGalL+k*m_nResGalL*m_nResGalKsi]<<endl;


	out.close();
}













////////////////////////////////////////////////////////////////////
//Density

CDensity::CDensity(){}
CDensity::~CDensity(){}

//Virutal methods of the CDensity class:
void CDensity::Norm(){}
void CDensity::Densify(){}
void CDensity::Zero(){}
RealType CDensity::GetDensitySpace(const SProbabilitySpacePointSpace SpacePoint){}
RealType CDensity::GetDensityRadio(const SProbabilitySpacePointRadio SpacePoint){}

RealType CDensity::GetDensity(const SProbabilitySpacePointRadio SpacePointRadio, const SProbabilitySpacePointSpace SpacePointSpace){}

void CDensity::AddPoint(const CPulsar Pulsar){}
void CDensity::AddValueAtPoint(const SProbabilitySpacePointRadio SpacePointRadio, const RealType fValue){}


void CDensity::Write(const char* FilePrefix){}
void CDensity::Check(){}
RealType CDensity::GetEfficiency() {}

RealType CDensity2x3D::GetEfficiency() 
{
	return (static_cast<RealType>(m_pSpace->m_nAddedPSRs) / m_pSpace->m_nTotalPSRs + static_cast<RealType>(m_pSpace->m_nAddedPSRs) / m_pSpace->m_nTotalPSRs)/2.;


}

CDensity2x3D::CDensity2x3D()
{
	m_pCfgPtr = CConfiguration::GetInstance();
	m_pRadio = new CRadioDensity(m_pCfgPtr->Plot.fZeroLogP, m_pCfgPtr->Plot.fZeroLogPdot, m_pCfgPtr->Plot.fZeroLogS1400, m_pCfgPtr->Plot.fMaxLogP, m_pCfgPtr->Plot.fMaxLogPdot, m_pCfgPtr->Plot.fMaxLogS1400, m_pCfgPtr->Plot.nResLogP, m_pCfgPtr->Plot.nResLogPdot, m_pCfgPtr->Plot.nResLogS1400 );
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		m_pSpace = new CSpaceDensity(m_pCfgPtr->Plot.fMaxDM, m_pCfgPtr->Plot.nResGalL,  m_pCfgPtr->Plot.nResGalKsi,  m_pCfgPtr->Plot.nResDM);
	}

}


CDensity2x3D::~CDensity2x3D()
{
	delete m_pRadio;
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		delete m_pSpace;
	}

}


void CDensity2x3D::AddValueAtPoint(const SProbabilitySpacePointRadio SpacePointRadio, const RealType fValue)
{
    //cout<<SpacePointRadio.nLogP<<" "<<SpacePointRadio.nLogPdot<<" "<<SpacePointRadio.nLogS1400<<" "<<fValue<<endl;
	m_pRadio->AddValueAtPoint(SpacePointRadio.nLogP,SpacePointRadio.nLogPdot,SpacePointRadio.nLogS1400,fValue);
}

RealType CDensity2x3D::GetDensityRadio(const SProbabilitySpacePointRadio SpacePoint)
{
	return m_pRadio->GetDensity(SpacePoint.nLogP,SpacePoint.nLogPdot,SpacePoint.nLogS1400);
}


RealType CDensity2x3D::GetDensitySpace(const SProbabilitySpacePointSpace SpacePoint)
{
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		return m_pSpace->GetDensity(SpacePoint.nGalL,SpacePoint.nGalKsi,SpacePoint.nDM);
	}
	else
	{
		return 1.;
	}
}

void CDensity2x3D::Check()
{
	cout<<m_pRadio->Check();
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		cout<<" "<<m_pSpace->Check();
	}
	cout<<endl;

}




RealType CDensity2x3D::GetDensity(const SProbabilitySpacePointRadio SpacePointRadio, const SProbabilitySpacePointSpace SpacePointSpace)
{
	return GetDensityRadio(SpacePointRadio)*GetDensitySpace(SpacePointSpace);
}

void CDensity2x3D::AddPoint(const CPulsar Pulsar)
{
	m_pRadio->AddPoint(Pulsar);
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		m_pSpace->AddPoint(Pulsar);
	}
}


void CDensity2x3D::Norm()
{
	m_pRadio->Norm();
	
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		m_pSpace->Norm();
	}
}


void CDensity2x3D::Densify()
{
	m_pRadio->Densify();
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		m_pSpace->Densify();
	}
}

void CDensity2x3D::Write(const char* FilePrefix)
{
	string Space(FilePrefix);
	Space+=".SpatialDensity.dat";
	string Radio(FilePrefix);
	Radio+=".RadioDensity.dat";
	
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		m_pSpace->Write(Space.c_str());
	}
	m_pRadio->Write(Radio.c_str());

}

void CDensity2x3D::Zero()
{
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D)
	{
		m_pSpace->Zero();
	}
	m_pRadio->Zero();

}


////////////////////////////////////////////
//Sparse

RealType CDensity2x3DSparse::GetEfficiency() 
{
	return (static_cast<RealType>(m_nAddedPSRsRadio) / m_nTotalPSRsRadio + static_cast<RealType>(m_nAddedPSRsSpace) / m_nTotalPSRsSpace)/2.;
}



CDensity2x3DSparse::CDensity2x3DSparse(const CObservations* pObservations)
{
	m_pCfgPtr=CConfiguration::GetInstance();
	CopyMapKeys(pObservations);

	m_fVolumeElementRadio=m_ComputeVolumeElementRadio();
	m_fVolumeElementSpace=m_ComputeVolumeElementSpace();
	Zero();
}

void CDensity2x3DSparse::CopyMapKeys(const CObservations* pObservations)
{

	for(std::map<SProbabilitySpacePointRadio, RealType>::const_iterator it = pObservations->MapRadio.begin(); it != pObservations->MapRadio.end(); ++it)
	{
		//efektywniejsze wkładnie - i tak są posortowane
		std::map<SProbabilitySpacePointRadio, RealType>::iterator endit =MapRadio.end();
		MapRadio.insert(endit,pair<SProbabilitySpacePointRadio,RealType>(it->first,0.));
	}

	for(std::map<SProbabilitySpacePointSpace, RealType>::const_iterator it = pObservations->MapSpace.begin(); it != pObservations->MapSpace.end(); ++it)
	{
		//efektywniejsze wkładnie - i tak są posortowane
		std::map<SProbabilitySpacePointSpace, RealType>::iterator endit =MapSpace.end();
		MapSpace.insert(endit,pair<SProbabilitySpacePointSpace,RealType>(it->first,0.));
	}

}

//	void m_AddToMap(const SProbabilitySpacePointRadio SpacePoint);
//
//jak jest element objetosci??
RealType CDensity2x3DSparse::m_ComputeVolumeElementSpace() const
{
	RealType dLogP( (m_pCfgPtr->Plot.fMaxLogP-m_pCfgPtr->Plot.fZeroLogP)/m_pCfgPtr->Plot.nResLogP );
	RealType dLogPdot( (m_pCfgPtr->Plot.fMaxLogPdot-m_pCfgPtr->Plot.fZeroLogPdot)/m_pCfgPtr->Plot.nResLogPdot);
	RealType dS1400( (m_pCfgPtr->Plot.fMaxLogS1400-m_pCfgPtr->Plot.fZeroLogS1400)/m_pCfgPtr->Plot.nResLogS1400 );
	return dLogP*dLogPdot*dS1400;
}

RealType CDensity2x3DSparse::m_ComputeVolumeElementRadio() const
{
	RealType dKsi( 2./m_pCfgPtr->Plot.nResGalKsi );
	RealType dL( 360./m_pCfgPtr->Plot.nResGalL );
	RealType dDM( m_pCfgPtr->Plot.fMaxDM/m_pCfgPtr->Plot.nResDM );
	return dKsi*dL*dDM;
}

void CDensity2x3DSparse::m_AddToMap(const SProbabilitySpacePointRadio Point)
{
	map<SProbabilitySpacePointRadio,RealType>::iterator RadioIterator;
	RadioIterator = MapRadio.find(Point);

	if( RadioIterator != MapRadio.end() )
	{
		RadioIterator->second += 1.; 
		m_nAddedPSRsRadio++;
	}
	m_nTotalPSRsRadio++;
}

void CDensity2x3DSparse::m_AddToMap(const SProbabilitySpacePointSpace Point)
{
	map<SProbabilitySpacePointSpace,RealType>::iterator SpaceIterator;
	SpaceIterator = MapSpace.find(Point);

	if( SpaceIterator != MapSpace.end() )
	{
		SpaceIterator->second += 1.;
		m_nAddedPSRsSpace++; 
	}
	m_nTotalPSRsSpace++;
}



void CDensity2x3DSparse::Zero()
{
	m_nAddedPSRsRadio=0;
	m_nTotalPSRsRadio=0;
	m_nAddedPSRsSpace=0;
	m_nTotalPSRsSpace=0;
	fDensityElement=1.;

	for(std::map<SProbabilitySpacePointSpace, RealType>::iterator it = MapSpace.begin(); it != MapSpace.end(); ++it)
	{
		it->second = 0.;
	}

	for(std::map<SProbabilitySpacePointRadio, RealType>::iterator it = MapRadio.begin(); it != MapRadio.end(); ++it)
	{
		it->second = 0.;
	}
}

RealType CDensity2x3DSparse::GetDensitySpace(const SProbabilitySpacePointSpace SpacePoint)
{
	return MapSpace[SpacePoint] * fDensityElement;

}

RealType CDensity2x3DSparse::GetDensityRadio(const SProbabilitySpacePointRadio SpacePoint)
{
	return MapRadio[SpacePoint] * fDensityElement;
}

RealType CDensity2x3DSparse::GetDensity(const SProbabilitySpacePointRadio SpacePointRadio, const SProbabilitySpacePointSpace SpacePointSpace)
{
	if( m_pCfgPtr->Sim.ModelDensityContainer == Density1x3DSparse )
	{
		return GetDensityRadio(SpacePointRadio);
	}
	else
	{
		return GetDensitySpace(SpacePointSpace)*GetDensityRadio(SpacePointRadio);
	}
}
	
void CDensity2x3DSparse::AddPoint(const CPulsar Pulsar)
{
	SProbabilitySpacePointRadio Radio = ComputeSpacePointCoordiantesRadio(Pulsar);
	m_AddToMap(Radio);
	
	SProbabilitySpacePointSpace Space = ComputeSpacePointCoordiantesSpace(Pulsar); 
	m_AddToMap(Space);
}


//Empty functions:
void CDensity2x3DSparse::Check(){}
void CDensity2x3DSparse::Densify()
{
	if( m_pCfgPtr->Sim.ModelDensityContainer == Density1x3DSparse )
	{
		fDensityElement /= m_fVolumeElementRadio;
	}
	else
	{
		fDensityElement /= m_fVolumeElementRadio*m_fVolumeElementSpace;

	}
}
void CDensity2x3DSparse::Norm()
{
    WARNING("Not adopted to the continous values");
	if( m_pCfgPtr->Sim.ModelDensityContainer == Density1x3DSparse )
	{
		fDensityElement /=m_nTotalPSRsRadio;
		
	}
	else
	{
		fDensityElement /=m_nTotalPSRsRadio*m_nTotalPSRsSpace;

	}


	//TODO błędne obliczanie normalizacji!!!
	//fDensityElement /= m_nTotalPSRsRadio;
}
void CDensity2x3DSparse::Write(const char* FilePrefix){}


///////////////////////////
//density4DSparse
RealType CDensity4DSparse::GetEfficiency() 
{
	return (static_cast<RealType>(m_nAddedPSRs) / m_nTotalPSRs );
}

RealType CDensity4DSparse::m_ComputeVolumeElement() const
{
	RealType dLogP( (m_pCfgPtr->Plot.fMaxLogP-m_pCfgPtr->Plot.fZeroLogP)/m_pCfgPtr->Plot.nResLogP );
	RealType dLogPdot( (m_pCfgPtr->Plot.fMaxLogPdot-m_pCfgPtr->Plot.fZeroLogPdot)/m_pCfgPtr->Plot.nResLogPdot);
	RealType dS1400( (m_pCfgPtr->Plot.fMaxLogS1400-m_pCfgPtr->Plot.fZeroLogS1400)/m_pCfgPtr->Plot.nResLogS1400 );
	RealType dDM( m_pCfgPtr->Plot.fMaxDM/m_pCfgPtr->Plot.nResDM );
	return dLogP*dLogPdot*dS1400*dDM;
}

void CDensity4DSparse::m_AddToMap(const SProbabilitySpacePoint4D Point)
{
	map<SProbabilitySpacePoint4D,RealType>::iterator Iterator;
	Iterator = Map4D.find(Point);

	if( Iterator != Map4D.end() )
	{
		Iterator->second += 1.; 
		m_nAddedPSRs++;
	}
	m_nTotalPSRs++;
}

void CDensity4DSparse::Zero()
{
	m_nAddedPSRs=0;
	m_nTotalPSRs=0;
	fDensityElement=1.;


	for(std::map<SProbabilitySpacePoint4D, RealType>::iterator it = Map4D.begin(); it != Map4D.end(); ++it)
	{
		it->second = 0.;
	}
}

void CDensity4DSparse::CopyMapKeys(const CObservations* pObservations)
{
	for(std::map<SProbabilitySpacePoint4D, RealType>::const_iterator it = pObservations->Map4D.begin(); it != pObservations->Map4D.end(); ++it)
	{
		//efektywniejsze wkładnie - i tak są posortowane
		std::map<SProbabilitySpacePoint4D, RealType>::iterator endit =Map4D.end();
		Map4D.insert(endit,pair<SProbabilitySpacePoint4D,RealType>(it->first,0.));
	}
}

RealType CDensity4DSparse::GetDensity(const SProbabilitySpacePoint4D SpacePoint)
{
	return Map4D[SpacePoint] * fDensityElement;
}

void CDensity4DSparse::AddPoint(const CPulsar Pulsar)
{
	SProbabilitySpacePoint4D Point = ComputeSpacePointCoordiantes4D(Pulsar);
	m_AddToMap(Point);
}
	
CDensity4DSparse::CDensity4DSparse(const CObservations* pObservations)
{
	m_pCfgPtr=CConfiguration::GetInstance();
	CopyMapKeys(pObservations);

	m_fVolumeElement=m_ComputeVolumeElement();
	Zero();
}
CDensity4DSparse::~CDensity4DSparse()
{
}

	//Empty functions
void CDensity4DSparse::Check()
{

}
void CDensity4DSparse::Norm()
{
    WARNING("Not adopted to the continous values");
	fDensityElement /= m_nTotalPSRs;
}
void CDensity4DSparse::Densify()
{
	fDensityElement /= m_fVolumeElement;
}
void CDensity4DSparse::Write(const char* FilePrefix)
{
}











void CDensityVector::Cumulate()
{
	pPeriodLogCD[0] = pPeriodLog[0] + fLogPLeft;
	pPeriodDotLogCD[0] = pPeriodDotLog[0] + fLogPdotLeft;
	pS1400LogCD[0] = pS1400Log[0] + fLogS1400Left;

	for(int i(1);i<nPeriodSize;++i) pPeriodLogCD[i]+=pPeriodLogCD[i-1]+pPeriodLog[i];
	for(int i(1);i<nPeriodDotSize;++i) pPeriodDotLogCD[i]+=pPeriodDotLogCD[i-1]+pPeriodDotLog[i];
	for(int i(1);i<nS1400Size;++i) pS1400LogCD[i] += pS1400LogCD[i-1]+pS1400Log[i];
}






void CDensityVector::AddPoint(const CPulsar PSR)
{
	nTotalNum++;

	AddPeriod(PSR.fPeriod);
	AddPeriodDot(PSR.fPeriodDot);
	AddS1400(PSR.fS1400);
}

void CDensityVector::AddPoint(const SCatalogueEntry OBS)
{
	nTotalNum++;

	AddPeriod(OBS.fP);
	AddPeriodDot(OBS.fPdot);
	AddS1400(OBS.fRadioFlux1400);
}


void CDensityVector::Zero()
{
	nTotalNum=0;
	fPeriodNumFraction=0.;
	fPeriodDotNumFraction=0.;
	fS1400NumFraction=0.;

	fLogPLeft=0.;
	fLogPRight=0.;
	fLogPdotLeft=0.;
	fLogPdotRight=0.;
	fLogS1400Left=0.;
	fLogS1400Right=0.;

	for(int i(0);i<nPeriodSize;++i) pPeriodLog[i]=0.;
	for(int i(0);i<nPeriodDotSize;++i) pPeriodDotLog[i]=0.;
	for(int i(0);i<nS1400Size;++i) pS1400Log[i]=0.;
	
	for(int i(0);i<nPeriodSize;++i) pPeriodLogCD[i]=0.;
	for(int i(0);i<nPeriodDotSize;++i) pPeriodDotLogCD[i]=0.;
	for(int i(0);i<nS1400Size;++i) pS1400LogCD[i]=0.;
}

void CDensityVector::Norm(const int nNum)
{
	if(nTotalNum>0)
	{
		fPeriodNumFraction/=static_cast<RealType>(nNum);
		fPeriodDotNumFraction/=static_cast<RealType>(nNum);
		fS1400NumFraction/=static_cast<RealType>(nNum);

		for(int i(0);i<nPeriodSize;++i) pPeriodLog[i]/=static_cast<RealType>(nNum);
		for(int i(0);i<nPeriodDotSize;++i) pPeriodDotLog[i]/=static_cast<RealType>(nNum);
		for(int i(0);i<nS1400Size;++i) pS1400Log[i]/=static_cast<RealType>(nNum);

		fLogPLeft/=static_cast<RealType>(nNum);
		fLogPRight/=static_cast<RealType>(nNum);
		fLogPdotLeft/=static_cast<RealType>(nNum);
		fLogPdotRight/=static_cast<RealType>(nNum);
		fLogS1400Left/=static_cast<RealType>(nNum);
		fLogS1400Right/=static_cast<RealType>(nNum);
	}
	Cumulate();
}

void CDensityVector::Norm()
{
//TODO add flag for Poisson
//    WARNING("Not adopted to the continous values");
//    PrintStacktrace();
	Norm(nTotalNum);
}

void CDensityVector::AddPeriod(const RealType Period)
{
	int nIndex = static_cast<int>((log10(Period)-m_pCfgPtr->Plot.fZeroLogP)/fDLogP);
	if(nIndex>=0 && nIndex<m_pCfgPtr->Plot.nResLogP) 
	{
		pPeriodLog[nIndex]++; 
		fPeriodNumFraction+=1.;
	}
	else
	{
		if(nIndex<0)
		{
			fLogPLeft++;
		}
		else
		{
			fLogPRight++;
		}

	}


}
void CDensityVector::AddPeriodDot(const RealType PeriodDot)
{
	int nIndex = static_cast<int>((log10(PeriodDot)-m_pCfgPtr->Plot.fZeroLogPdot)/fDLogPdot);
	if(nIndex>=0 && nIndex<m_pCfgPtr->Plot.nResLogPdot)
	{
		pPeriodDotLog[nIndex]++; 
		fPeriodDotNumFraction+=1.;
	}
	else
	{
		if(nIndex<0)
		{
			fLogPdotLeft++;
		}
		else
		{
			fLogPdotRight++;
		}

	}


}
void CDensityVector::AddS1400(const RealType S1400)
{
	int nIndex = static_cast<int>((log10(S1400)-m_pCfgPtr->Plot.fZeroLogS1400)/fDS1400);
	if(nIndex>=0 && nIndex<m_pCfgPtr->Plot.nResLogS1400) 
	{
		pS1400Log[nIndex]++; 
		fS1400NumFraction+=1.;
	}
	else
	{
		if(nIndex<0)
		{
			fLogS1400Left++;
		}
		else
		{
			fLogS1400Right++;
		}
	}


}

CDensityVector::CDensityVector()
{
	m_pCfgPtr = CConfiguration::GetInstance();

	fDLogP = (m_pCfgPtr->Plot.fMaxLogP-m_pCfgPtr->Plot.fZeroLogP)/m_pCfgPtr->Plot.nResLogP ;
	fDLogPdot = (m_pCfgPtr->Plot.fMaxLogPdot-m_pCfgPtr->Plot.fZeroLogPdot)/m_pCfgPtr->Plot.nResLogPdot;
	fDS1400 = (m_pCfgPtr->Plot.fMaxLogS1400-m_pCfgPtr->Plot.fZeroLogS1400)/m_pCfgPtr->Plot.nResLogS1400 ;

	nPeriodSize = m_pCfgPtr->Plot.nResLogP;
	nPeriodDotSize = m_pCfgPtr->Plot.nResLogPdot;
	nS1400Size = m_pCfgPtr->Plot.nResLogS1400;



	pPeriodLog = new RealType[m_pCfgPtr->Plot.nResLogP];
	pPeriodDotLog = new RealType[m_pCfgPtr->Plot.nResLogPdot];
	pS1400Log = new RealType[m_pCfgPtr->Plot.nResLogS1400];
	
	pPeriodLogCD = new RealType[m_pCfgPtr->Plot.nResLogP];
	pPeriodDotLogCD = new RealType[m_pCfgPtr->Plot.nResLogPdot];
	pS1400LogCD = new RealType[m_pCfgPtr->Plot.nResLogS1400];

	Zero();
}

CDensityVector::~CDensityVector()
{
	delete [] pPeriodLog;
	delete [] pPeriodDotLog;
	delete [] pS1400Log;
	
	delete [] pPeriodLogCD;
	delete [] pPeriodDotLogCD;
	delete [] pS1400LogCD;
}



CDensity3x1D::CDensity3x1D()
{
	m_pCfgPtr = CConfiguration::GetInstance();
	Zero();
}

CDensity3x1D::~CDensity3x1D()
{
}

void CDensity3x1D::Zero()
{
	Vector.Zero();
}

void CDensity3x1D::AddPoint(const CPulsar Pulsar)
{
	Vector.AddPoint(Pulsar);
}

void CDensity3x1D::Densify()
{

}

void CDensity3x1D::Norm()
{
	Vector.Norm();
}



RealType DValue3D(const CDensityVector *pVA, const CDensityVector *pVB)
{
	RealType fDPLog(0.);
	RealType fDPDotLog(0.);
	RealType fDS1400Log(0.);

	CConfiguration *pCfgPtr = CConfiguration::GetInstance();

	for(int i(0);i< pCfgPtr->Plot.nResLogP  ;i++) 
	{
		fDPLog=fmax(fDPLog,fabs(pVA->pPeriodLogCD[i] - pVB->pPeriodLogCD[i]));
		//cout<<i<<" "<<pVA->pPeriodLog[i]<<" "<<pVB->pPeriodLog[i]<<endl;
	}

	for(int i(0);i< pCfgPtr->Plot.nResLogPdot;i++) 
	{
		fDPDotLog=fmax(fDPDotLog,fabs(pVA->pPeriodDotLogCD[i] - pVB->pPeriodDotLogCD[i]));
	}
	for(int i(0);i< pCfgPtr->Plot.nResLogS1400;i++) 
	{
		fDS1400Log=fmax(fDS1400Log,fabs(pVA->pS1400LogCD[i] - pVB->pS1400LogCD[i]));
	}
	//RealType fPeriodNorm = pVA->fPeriodNumFraction * pVB->fPeriodNumFraction;
	//RealType fPeriodDotNorm = pVA->fPeriodDotNumFraction * pVB->fPeriodDotNumFraction;
	//RealType fS1400Norm = pVA->fS1400NumFraction * pVB->fS1400NumFraction;

	return sqrt(fDPLog*fDPLog + fDPDotLog*fDPDotLog + fDS1400Log*fDS1400Log);
}





