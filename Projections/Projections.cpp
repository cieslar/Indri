
#include<iostream>
#include<fstream>
#include<Statistics.hpp>
#include<cmath>
#include<Debug.hpp>
#include<cstring>
#include<Randomizer.hpp>

#include<Projections.hpp>
#include<iomanip>

using namespace std;




bool DoubleEq(const double fLeft, const double fRight)
{
    double fRelativeEps = fLeft * std::numeric_limits<double>::epsilon();
    if( fabs(fLeft-fRight) <= fabs(fRelativeEps)) return true;//FIXME czy to jest prawda?
    return false;
    
}





double& CDoubleArray2D::Data(const int nX, const int nY)
{
    if(!(nX<nWidth && nY<nHight))
    {
	ERROR("Memmory access out of range");
    }
    return pData[nX+nY*nWidth];
}


CDoubleArray2D::CDoubleArray2D(const int nArrayWidth, const int nArrayHight)
{
    nWidth=nArrayWidth;
    nHight=nArrayHight;
    nNum=nWidth*nHight;
    nMemMax=ComputeMem();

    pData=new double[nNum];
    Zero();
}

CDoubleArray2D::~CDoubleArray2D()
{
	delete [] pData;
}

void CDoubleArray2D::BinWrite(const char* FileName) const
{
	std::ofstream out(FileName, std::ios::binary);
	out.write((char *)pData,nMemMax);
	out.close();
}
void CDoubleArray2D::BinRead(const char* FileName)
{
	std::ifstream in(FileName, std::ios::binary);
	in.read((char *)pData,nMemMax);
	in.close();
}

void CDoubleArray2D::Write(const char* FileName) const
{
	std::ofstream out(FileName);
	for(int i(0);i<nWidth;i++)
	{
	    for(int j(0);j<nHight;j++)
	    {
		out<<pData[i+j*nWidth]<<" ";
	    }
	    out<<std::endl;
	}
	out.close();
}
void CDoubleArray2D::Read(const char* FileName)
{
	std::ifstream in(FileName);
	for(int i(0);i<nWidth;i++)
	{
	    for(int j(0);j<nHight;j++)
	    {
		in>>Data(i,j);
	    }
	}
	in.close();
}

void CDoubleArray2D::Zero()
{
	double *p = pData, *last = pData + nNum;
	while(p != last) *(p++) = 0.;


	//memset(pData,0.,nMemMax);
}

size_t  CDoubleArray2D::ComputeMem()
{
	return nNum*sizeof(double);
}















/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


CProjection::CProjection(const int nResX, const double fLZero) : 
    m_SquareProjection(nResX,nResX/2)
{
	if((nResX/2)*2 != nResX)
	{
		WARNING("Resolution is not even.");
	}
	m_nResX=nResX;
	m_nResY=m_nResX/2;
	m_fZeroX=M_PI;
	m_fZeroY=M_PI/2.;

	m_fDeltaX = 2.*M_PI/m_nResX;
	m_fDeltaY =    M_PI/m_nResY;

}

void CProjection::AddPoint(const double fL, const double fB)
{
	if((fabs(fB)>M_PI/2.) || (fabs(fL)>M_PI))
	{
		ERROR("Galactic coodrinates out of the domain");
	}
	m_SquareProjection.Data(static_cast<int>((fL+m_fZeroX)/m_fDeltaX),static_cast<int>((fB+m_fZeroY)/m_fDeltaY))++;
}

void CProjection::SphereDensify()
{
	double fFi(0.),fTheta(0.),fRadius(m_nResX/2.);
	for(int i(0);i<m_nResX;i++)
	{
		fFi=i*m_fDeltaX - m_fZeroX;
		for(int j(0);j<m_nResY;j++)
		{
			fTheta=j*m_fDeltaY - m_fZeroY;
			m_SquareProjection.Data(i,j)/= fRadius*fRadius * cos(fTheta) * m_fDeltaX * m_fDeltaY;
		}
	}
}

void CProjection::Norm()
{
	double fTotal(0.);
	for(int i(0);i<m_nResX;i++)
	{
		for(int j(0);j<m_nResY;j++)
		{
			fTotal+=m_SquareProjection.Data(i,j);
		}
	}

	for(int i(0);i<m_nResX;i++)
	{
		for(int j(0);j<m_nResY;j++)
		{
			m_SquareProjection.Data(i,j)/=fTotal;
		}
	}
}



void CProjection::Zero()
{
	m_SquareProjection.Zero();

}

void CProjection::Write(const char *FileName)
{
	std::ofstream out(FileName);
	for(int j(0);j<m_nResY;j++)
	{
		for(int i(0);i<m_nResX;i++)
		{
			out<<m_SquareProjection.Data(i,j)<<" ";
		}
		out<<std::endl;
	}
	out.close();

}


void CProjection::BinWrite(const char *FileName)
{
	std::ofstream out(FileName, std::ios::binary);
	out.write((char *)m_SquareProjection.pData,m_SquareProjection.nMemMax);
	out.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




double ThetaFunc(const double fTheta, const double fB)
{
	return (2.*fTheta + sin(2.*fTheta) - M_PI * sin(fB));
}




double CMollweideProjection::m_ComputeTheta(const double fB)
{
    	//FIXME zamienic na porownanie doublowe
    	if(DoubleEq(fB, M_PI/2.)) return M_PI/2.;
    	if(DoubleEq(fB,-M_PI/2.)) return -M_PI/2.;
    	double fThetaP(0.), fThetaN(0.);
    	
	fThetaN=M_PI*sin(fB)/2.;
	double fEps(std::numeric_limits<double>::epsilon()), fRelativeEps,fDelta;
  	do
    	{
		fThetaP = fThetaN;
		fThetaN = fThetaP - (2.*fThetaP + sin(2.*fThetaP) - M_PI*sin(fB)) / 16.*(2. + 2.*cos(2.*fThetaP));
    		//fRelativeEps = fThetaN *100.* fEps;
		fDelta=fabs(fThetaN-fThetaP);
    	}
	    while(!(fDelta<=1e-6));
    
	return fThetaN;

}


void CMollweideProjection::Project()
{		
	double fL(0.), fB(0.), fTheta(0.), fX(0.), fY(0.);
	const double sq2 = sqrt(2.);
	double fLZero(0.);

	m_fRadius=1.;
	for(int i(0);i<m_nResX;i++)
	{
		fL=i*m_fDeltaX - m_fZeroX;
		for(int j(0);j<m_nResY;j++)
		{
			fB=j*m_fDeltaY - m_fZeroY;
			fTheta = m_ComputeTheta(fB);
			fX = m_fRadius * (2.*sq2/M_PI) * (fL-fLZero) * cos(fTheta);
			fY = m_fRadius * sq2 * sin(fTheta);

			//cout<<fX<<" "<<fY<<" "<<static_cast<int>((fL+m_fZeroX)/m_fDeltaX)<<" "<<static_cast<int>((fB+m_fZeroY)/m_fDeltaY)<<endl;
			//m_MollweideProjection.Data(static_cast<int>((fL+m_fZeroX)/m_fDeltaX),static_cast<int>((fB+m_fZeroY)/m_fDeltaY)) += m_SquareProjection.Data(i,j);
			m_MollweideProjection.Data(static_cast<int>((fX+m_fZeroX)/m_fDeltaX),static_cast<int>((fY+m_fZeroY)/m_fDeltaY)) += 10.;
		}
	}
}


#if 0
void CMollweideProjection::Project()
{		
	double fL(0.), fB(0.), fTheta(0.), fX(0.), fY(0.);
	const double sq2 = sqrt(2.);
	double fLZero(0.);

	m_fRadius=1.;

	for(double a(-M_PI);a<=M_PI;a+=1e-1) for(double b(-M_PI/2.);b<=M_PI/2.;b+=1e-1)	
	{
		for(double b(-M_PI/2.);b<=M_PI/2.;b+=1e-1)	
		{
			fTheta = m_ComputeTheta(b);
			fX = m_fRadius * (2.*sq2/M_PI) * (a) * cos(fTheta);
			fY = m_fRadius * sq2 * sin(fTheta);


			cout<<fX<<" "<<fY<<endl;
		}
		cout<<endl;
	}
}
#endif

void CMollweideProjection::Write(const char *FileName)
{
	std::ofstream out(FileName);
	for(int j(0);j<m_nResY;j++)
	{
		for(int i(0);i<m_nResX;i++)
		{
			out<<m_MollweideProjection.Data(i,j)<<" ";
		}
		out<<std::endl;
	}
	out.close();

}


void CMollweideProjection::BinWrite(const char *FileName)
{
	std::ofstream out(FileName, std::ios::binary);
	out.write((char *)m_MollweideProjection.pData,m_SquareProjection.nMemMax);
	out.close();
}


void CMollweideProjection::Zero()
{
	m_MollweideProjection.Zero();
	m_SquareProjection.Zero();
}



#if 0

void MollweideTest()
{
	CRandomizer Ran;
	double fL(0.), fB(0.);

	CMollweideProjection Proj(2048);
	for(int i(0); i<10000000; i++)
	{
//		fL=Ran.Flat()*2.*M_PI - M_PI;
//		fB=Ran.Flat()*M_PI - M_PI/2.;
//		Proj.AddPoint(fL,fB);
	}
	Proj.SphereDensify();
	Proj.Norm();
	Proj.Project();
	Proj.Write("MollweideProj.dat");


}


int main()
{
	MollweideTest();
	return 0;
}


#endif
