#include<Randomizer.hpp>
#include<Debug.hpp>

#include<cmath>
#include<cstdlib>
#include<time.h>

#include<iostream>


using namespace std;

unsigned int CRandomizer::RandomInt(const unsigned int nModulo)
{
	m_nCalls++;
	//FIXME czy musze castowac??
	return static_cast<unsigned int>(m_pMTR->randInt()) % nModulo;
}

CRandomizer::~CRandomizer()
{
//	cout<<"RanCalls: "<<nRanCount<<endl<<"MathCalls: "<<nMathCount<<endl;
	delete m_pMTR;

	if(m_pfGaussCDF!=NULL) delete [] m_pfGaussCDF;
}

CRandomizer::CRandomizer()
{
//	nRanCount=nMathCount=nTotCalls=0;
	m_bGaussCDF=false;

//	m_nCallsInt=0;
	m_nCalls=0;
	Seed(static_cast<unsigned int>(GetTimestamp()%1000000000));
//	m_Seed(time(NULL));
}

RealType CRandomizer::Flat()
{
	return Random();
}


RealType CRandomizer::FlatOpen()
{
	RealType a;
	do
	{
		a=Random();
		//cout<<a<<endl;
	}
	while(a==0. || a==1.);
	return a;
}



RealType CRandomizer::FlatLeftOpen()
{
	RealType a;
	do
	{
		a=Random();
		//cout<<a<<endl;
	}
	while(a==0.);
	return a;
}
RealType CRandomizer::FlatMinMax(RealType fMin, RealType fMax)
{
	return Flat()*(fMax-fMin) + fMin;
}

RealType CRandomizer::LogFlatMinMax(RealType fMin, RealType fMax)
{
	return pow(10., Flat()*(log10(fMax)-log10(fMin)) + log10(fMin) );
}

RealType CRandomizer::FlatMeanSigma(RealType fMean, RealType fSigma)
{
	return (Flat()*2. -1.) * ( sqrt(3.) * fSigma ) + fMean;

}

RealType CRandomizer::LogFlatMeanSigma(RealType fMean, RealType fSigma)
{
	return pow(10., (Flat()*2. -1.) * ( sqrt(3.) * log10(fSigma) ) + log10(fMean));
}

RealType CRandomizer::LogFlatMeanLogSigma(RealType fMean, RealType fLogSigma)
{
	return pow(10., (Flat()*2. -1.) * ( sqrt(3.) * fLogSigma ) + log10(fMean));
}

RealType CRandomizer::PositiveFlatMeanSigma(RealType fMean, RealType fSigma)
{ 
	RealType fResult; 
	do
	{
		fResult = FlatMeanSigma(fMean, fSigma);
	}
	while(fResult<0.);
	return fResult;
}

RealType CRandomizer::PositiveLogFlatMeanSigma(RealType fMean, RealType fSigma)
{	
	RealType fResult;
	do
	{
		fResult = pow(10., (Flat()*2. -1.) * ( sqrt(3.) * log10(fSigma) ) + log10(fMean));
	}
	while(fResult<0.);
	return fResult;
}
RealType CRandomizer::PositiveLogFlatMeanLogSigma(RealType fMean, RealType fLogSigma)
{	
	RealType fResult;
	do
	{
		fResult = pow(10., (Flat()*2. -1.) * ( sqrt(3.) * fLogSigma ) + log10(fMean));
	}
	while(fResult<0.);
	return fResult;
}



void CRandomizer::Seed(const unsigned int nSeed, const unsigned int nRollCount)
{
	m_nSeed=nSeed;
	//dynamiczna tablica aby nie ciazyla
        //jeden seed -> okres 2^32
	//pelna tablica seedow -> okres 2^19937-1
        MTRand::uint32 *pSeed = new MTRand::uint32[MTRand::N];

        for( int i(0); i < MTRand::N; ++i ) pSeed[i] = m_nSeed * i;

	if(!m_pMTR) m_pMTR=new MTRand(pSeed);
	else m_pMTR->seed(pSeed);

	delete [] pSeed;
}

RealType CRandomizer::Random()
{
	m_nCalls++; 
#if 0
	RealType result;
	//FIXME podpiac pod singleton configowy:
	switch(m_pCfgPtr->nRandomizerType)
	{
		case MTRRandomizer:
			result=static_cast<RealType>(m_pMTR->rand());
			break;
		default:
			result=static_cast<RealType>(rand())/RAND_MAX;
			break;
	}
#endif
	return static_cast<RealType>(m_pMTR->rand());
}

void CRandomizer::PrintInfo()
{
	cout<<"Randomizer info:"<<endl;
	cout<<"type:  1"<<endl;
	//cout<<"type:  "<<m_pCfgPtr->nRandomizerType<<endl;
	cout<<"seed:  "<<m_nSeed<<endl;
	cout<<"calls: "<<m_nCalls<<endl;
}


RealType CRandomizer::m_GaussDistribution(const RealType fX, const RealType fSigma, const RealType fMi)
{
	return (1./ (fSigma*sqrt(2.*M_PI)) )*exp(-(fX-fMi)*(fX-fMi)/(2.*fSigma*fSigma));
}

#if 0
RealType CRandomizer::m_LogNormalDistribution(const RealType fX, const RealType fSigma, const RealType fMi)
{
	return 0.;
	//FIXME! to jest dosc mocno schrznieone fMi jest tej samej sednostki co fX!
	//return (1./ (fX*fSigma*sqrt(2.*M_PI)) )*exp(-(log(fX)-fMi)*(log(fX)-fMi)/(2.*fSigma*fSigma));
}
#endif

RealType CRandomizer::m_ExponentialDistribution(const RealType fX, const RealType fLambda)
{
	return fLambda*exp(-fLambda*fX);
}



RealType CRandomizer::m_BimodalGaussDistribution(const RealType fX, const RealType fWeight, const RealType fSigma1, const RealType fSigma2, const RealType fMi1, const RealType fMi2)
{
	return 	(fWeight/(sqrt(2.*M_PI)*fSigma1)) * exp(-(fX-fMi1)*(fX-fMi1)/(2.*fSigma1*fSigma1)) + \
		((1.-fWeight)/(sqrt(2.*M_PI)*fSigma2)) * exp(-(fX-fMi2)*(fX-fMi2)/(2.*fSigma2*fSigma2));
}

RealType CRandomizer::BimodalGauss(const RealType fWeight, const RealType fSigma1, const RealType fSigma2, const RealType fMi1, const RealType fMi2)
{
	RealType fMax, fX, fPhi;
	RealType fSigma= (fSigma1 > fSigma2) ? fSigma1 : fSigma2;
	fMax=m_BimodalGaussDistribution(0., fWeight, fSigma1, fSigma2, fMi1, fMi2);
	do
	{
		fX=Flat()*6.*fSigma-3.*fSigma;
		fPhi=Flat()*fMax;
	}
	while(fPhi>m_BimodalGaussDistribution(fX, fWeight, fSigma1, fSigma2, fMi1, fMi2));
	//cout<<fX<<" "<<fPhi<<endl;
	return fX;
}





//Draw random number from -3sigma + 3sigma from Normal Distribution
//Will be slow :(
RealType CRandomizer::GaussBruteForce(const RealType fSigma, const RealType fMi)
{
	RealType fMax, fX, fPhi;

	fMax=m_GaussDistribution(fMi, fSigma, fMi);
	do
	{
		fX=Flat()*6.*fSigma-3.*fSigma +fMi;
		fPhi=Flat()*fMax;
	}
	while(fPhi>m_GaussDistribution(fX,fSigma,fMi));
	//cout<<fX<<" "<<fPhi<<endl;
	return fX;
}
Cartesian CRandomizer::Gauss2D(const RealType fSigma, const RealType fMi)
{
{
	RealType fMax, fX, fY,fPhi;
	fMax=m_GaussDistribution(fMi, fSigma, fMi)*m_GaussDistribution(fMi, fSigma, fMi);
	do
	{
		fX=Flat()*6.*fSigma-3.*fSigma +fMi;
		fY=Flat()*6.*fSigma-3.*fSigma +fMi;
		fPhi=Flat()*fMax;
	}
	while(fPhi>m_GaussDistribution(fX,fSigma,fMi)*m_GaussDistribution(fY,fSigma,fMi));
	//cout<<fX<<" "<<fPhi<<endl;
	Cartesian Result={fX,fY,0.};
	return Result;
}





}



#if 0
RealType CRandomizer::LogNormal(const RealType fSigma, const RealType fMi)
{
	RealType fMax, fX, fPhi;

	fMax=m_LogNormalDistribution(fMi, fSigma, fMi);//Czy w zerze?
	do
	{
		do
		{
			fX=Flat()*6.*fSigma-3.*fSigma +fMi;//FIXME - ile sigm potrzebuje w lognorm?
		}
		while(fX<0.);

		fPhi=Flat()*fMax;
	}
	while(fPhi>m_LogNormalDistribution(fX,fSigma,fMi));
	//cout<<fX<<" "<<fPhi<<endl;
	return fX;
}
#endif



RealType CRandomizer::PositiveGaussBruteForce(const RealType fSigma, const RealType fMi)
{
	RealType fMax, fX, fPhi;

	fMax=m_GaussDistribution(fMi, fSigma, fMi);
	do
	{	
		do
		{
			fX=Flat()*6.*fSigma-3.*fSigma +fMi;
		}
		while(fX<0.);
		fPhi=Flat()*fMax;
	}
	while(fPhi>m_GaussDistribution(fX,fSigma,fMi));
	//cout<<fX<<" "<<fPhi<<endl;
	return fX;
}

RealType CRandomizer::Exponential(const RealType fLambda)
{
	RealType fMax, fX, fPhi;

	fMax=m_ExponentialDistribution(0., fLambda);
	do
	{
		fX=Flat()*3./fLambda;//FIXME Doczytac gdzie lezy 90% w wykladniczym
		fPhi=Flat()*fMax;
	}
	while(fPhi>m_ExponentialDistribution(fX,fLambda));
	//cout<<fX<<" "<<fPhi<<endl;
	return fX;
}

RealType CRandomizer::DoubleSidedExponential(const RealType fLambda)
{
	RealType fMax, fX, fPhi;

	fMax=0.5*m_ExponentialDistribution(0., fLambda);
	do
	{
		fX=Flat()*6./fLambda-3./fLambda;
		fPhi=Flat()*fMax;
	}
	while(fPhi>0.5*m_ExponentialDistribution(((fX < 0.) ? -fX : fX),fLambda));
	//cout<<fX<<" "<<fPhi<<endl;
	return fX;
}



#include <cstdlib>
#include <limits>
#include<fstream>
RealType CRandomizer::GenerateGaussianBoxMuller(const RealType fSigma, const RealType fMi)
{
	if (fSigma<0.) 
	{
		//ERROR("Sigma in Gauss distribution is < 0");
		WARNING("Sigma in Gauss distribution is < 0");
		PrintStacktrace();
	}
	const RealType epsilon = std::numeric_limits<RealType>::min();
	const RealType two_pi = 2.0*3.14159265358979323846;
 
	static RealType z0, z1;
//	TTimestamp B0=GetTimestamp();
//	static bool generate;
//	generate = !generate;
 
//	if (!generate)
//	   return z1 * sigma + mu;
 
	RealType u1, u2;
	do
	 {
	   u1 = Flat();
	   u2 = Flat();
	 }
	while ( u1 <= epsilon );
 
//	TTimestamp B1=GetTimestamp();
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
//	TTimestamp B2=GetTimestamp();
//	nRanCount += B1-B0;
//	nMathCount += B2-B1;
	//z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * fSigma + fMi;
}


RealType CRandomizer::Gauss(const RealType fSigma, const RealType fMi)
{
	return GenerateGaussianBoxMuller(fSigma,fMi);
}

RealType CRandomizer::PositiveGauss(const RealType fSigma, const RealType fMi)
{
	RealType fRes;
	do
	{
		fRes=GenerateGaussianBoxMuller(fSigma,fMi);
	}
	while(fRes<0.);
	return fRes;
}



void TestRandomizer()
{
	CRandomizer R1;
	RealType fRes;
	TTimestamp A0, A1, B0,B1;
	A0=GetTimestamp();
	int n(1000000);
	for(int i(0);i<n;i++) 
	{
		R1.GenerateGaussianBoxMuller(1.,0.);
	}
	A1=GetTimestamp();
	B0=GetTimestamp();
	for(int i(0);i<n;i++) 
	{
		R1.GaussBruteForce(1.,0.);
	}
	B1=GetTimestamp();
	cout<<A1-A0<<endl<<B1-B0<<endl;
}



RealType CRandomizer::InverseTransformSamplingGauss(const RealType fSigma, const RealType fMi)
{
	if(!m_bGaussCDF) m_PrepGaussCDF(m_fGaussCDFResolution);

	RealType fX = Flat()*m_fSigmaToleranceInterval;
	RealType fSign = (Flat() >0.5) ? -1. : 1.;

	int nRight(m_nGaussCDFSize-1);
	int nLeft(0);
	int nMiddle;

//	while(nLeft<m_nGaussCDFSize && fX >= m_pfGaussCDF[nLeft]) ++nLeft;

//Sprawdzic ktore realnie jest szybsze jeśli chodzi o oczekiwany przedzial dla fX.
//	RealType fK=(fX- m_pfGaussCDF[nLeft])/( m_pfGaussCDF[nLeft-1] - m_pfGaussCDF[nLeft]);

	while(nLeft<nRight)
	{
		if( (nRight-nLeft) <= 1) break;

		nMiddle = ( nLeft + nRight ) / 2;

		if( fX <= m_pfGaussCDF[nMiddle] )
		{
			nRight = nMiddle;
		}
		else
		{
			nLeft = nMiddle + 1;
		}
	}	
	double fK(0.);
	if(nLeft<nRight) fK=(fX- m_pfGaussCDF[nLeft])/( m_pfGaussCDF[nRight] - m_pfGaussCDF[nLeft]);
	return fSign*((nLeft+fK)*m_fGaussCDFStep)*fSigma + fMi;
}




void CRandomizer::m_PrepGaussCDF(const RealType fStep)
{
	m_fGaussCDFStep=fStep;
	m_nGaussCDFSize = static_cast<int>(m_nSigmaRange/m_fGaussCDFStep) + 2;
	m_pfGaussCDF = new RealType[m_nGaussCDFSize];
	double fFrom, fTo;

	//m_pfGaussCDF[0]=m_RombergIntegralGauss(0.,m_fGaussCDFStep )+0.5 ;
	m_pfGaussCDF[0]=0. ;
	for(int i(1);i<m_nGaussCDFSize;i++)
	{
		fFrom =   (i-1)   * m_fGaussCDFStep;
		fTo   = (i) * m_fGaussCDFStep;
		m_pfGaussCDF[i]=m_RombergIntegralGauss(fFrom, fTo) + m_pfGaussCDF[i-1];
	}


	RealType fTotal(m_pfGaussCDF[m_nGaussCDFSize-1]);
	for(int i(0);i<m_nGaussCDFSize;i++) m_pfGaussCDF[i]/=fTotal;


	m_bGaussCDF=true;
}



RealType CRandomizer::m_RombergIntegralGauss(const RealType fFrom, const RealType fTo)
{
	const int nSize(4);
	double Tab[nSize];

	//Starting number of intervals
	int nNumber(4);
	double fH( (fFrom-fTo)/nNumber );

	//First integral
	Tab[0] = (fH/2.) * ( m_GaussDistribution(fFrom) + m_GaussDistribution(fTo));

	for(int i(1); i<nNumber; i++)
	{
		Tab[0] += fH*m_GaussDistribution(fFrom + i*fH);
	}

	//Rest of integrals
	for(int i(1);i<nSize;i++)
	{	
		Tab[i] = Tab[i-1]/2.;
		for(int j(0); j<nNumber; j++)
		{
			Tab[i] += (fH/2.) * m_GaussDistribution(fFrom + fH/2. + j*fH);
		}

		fH/=2.;
		nNumber*=2;
	}

	//In-place Richardson Extrapolation
	double fK(1.);
	for(int i(1);i<nSize;i++)
	{
		fK*=4.;
		for(int j(0);j<nSize-i;j++)
		{
			Tab[j] = (1./(fK-1.))*(fK*Tab[j+1] -Tab[j]);
		}
	}

	//Return the best result
	return Tab[0];
}


bool CRandomizer::RandomLogical()
{
	if(Flat()>0.5) return true;
	return false;
}


RealType CRandomizer::m_MaxwellDistribution(const RealType fX, const RealType fA)
{
	return sqrt(2./M_PI) * (fX*fX * exp(- fX*fX / (2. * fA*fA) ) / (fA*fA*fA));
}

RealType CRandomizer::Maxwell(const RealType fA)
{
	const RealType fEps(1e-4);
	RealType fMaxY( m_MaxwellDistribution( sqrt(2.)*fA ,fA) + fEps);
	RealType fSigma(sqrt(fA*fA*(3.*M_PI-8.)/M_PI));
	RealType fX(0.),fY(0.);
	RealType fMaxX(sqrt(2.)*fA + 5.*fSigma);
	//cout<<fA<<" "<<fMaxY<<" "<<fSigma<<" "<<fMaxX<<endl;
	do
	{
		fX=Flat()*fMaxX;
		fY=Flat()*fMaxY;
	}
	while( fY > m_MaxwellDistribution(fX,fA) );

	return fX;
}


Cartesian CRandomizer::UniformShpereDirection()
{
	//wylosowanie kierunku
	//http://mathworld.wolfram.com/SpherePointPicking.html
	RealType x1,x2;
	do
	{
		x1=FlatOpen()*2. - 1.;
		x2=FlatOpen()*2. - 1.;

	}
	while(x1*x1 + x2*x2 >= 1.);

	///Unit vector for the direction
	Cartesian Cart;

	Cart.fX=2.*x1*sqrt(1.-x1*x1 -x2*x2);
	Cart.fY=2.*x2*sqrt(1.-x1*x1 -x2*x2);
	Cart.fZ=1.-2.*(x1*x1 + x2*x2);

	return Cart;
}


