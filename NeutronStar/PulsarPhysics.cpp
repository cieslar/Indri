#include<PulsarPhysics.hpp>

#include<ConfigContainer.hpp>

#include<cmath>
#include<limits>
#include<fstream>
using namespace std;

RealType MYears2Seconds(const RealType fMYears)
{
	return 1e6*365.*24.*60.*60.*fMYears;//FIXME zastanowic sie jaki rok przyjac ;)
}

RealType Seconds2MYears(const RealType fSeconds)
{
	return fSeconds/(1e6*365.*24.*60.*60.);
}



void CPulsarPhysics::Evolve(const RealType fTimeStep, const bool bOneStep)
{
#if 0
	if(bOneStep)
	{
		//cout<<fAge<<" "<<fTimeStep<<endl;
		int i(256);
		while (fTimeStep/static_cast<RealType>(i) > 10.) i<<=1;
	//while (fTimeStep/static_cast<RealType>(i) > m_pCfgPtr->Sim.fMaxTimeStep) i<<=1;
	//if(m_pCfgPtr->Sim.bVerbose) cout<<fTimeStep<<" "<<i<<endl;
		for (int j(0);j<i;j++) 
		{
			EvolveRK4Adaptive(fTimeStep/static_cast<RealType>(i));
	//EvolveRK4(fTimeStep);

			//fAge+=fTimeStep/static_cast<RealType>(i);
			//fRelativeAge+=fTimeStep/static_cast<RealType>(i);
			//m_ComputePPdotBField(fTimeStep);
		}
	}
#endif
	//EvolveRK4Adaptive(fTimeStep);
	m_ComputePPdotBField(fTimeStep);
}


void CPulsarPhysics::EvolveForwardNewton(const RealType fTimeStep)
{
	fPeriod += fPeriodDot * MYears2Seconds(fTimeStep);//FIXME na RK4
	m_EvolveBField(fTimeStep);
	fPeriodDot=m_PeriodDerivative(fPeriod,fBField);
}



void CPulsarPhysics::EvolveRK4(const RealType fTimeStep)
{
//BField nie jest w tym prostym modelu uwiklane z P i moge potraktowac je jako funkcje t. Inaczej potrzebna bedzie wersja wektorowa.
	RealType fTimeStepSeconds =  MYears2Seconds(fTimeStep);

	RealType fA = fPeriodDot;//wyliczone w poprzednim kroku lub w init'cie
	RealType fB = m_PeriodDerivative(fPeriod + fA*fTimeStepSeconds/2., m_NextBField(fBField, fTimeStep/2.));
	RealType fC = m_PeriodDerivative(fPeriod + fB*fTimeStepSeconds/2., m_NextBField(fBField, fTimeStep/2.));
	fBField = m_NextBField(fBField, fTimeStep);//i tak bede zapisywal B z konca kroku
	RealType fD = m_PeriodDerivative(fPeriod + fC*fTimeStepSeconds   , fBField);
	
	fPeriod += (fA + 2.*fB + 2.*fC + fD)*fTimeStepSeconds/6.;

	fPeriodDot=m_PeriodDerivative(fPeriod,fBField);
}


void CPulsarPhysics::EvolveRK4Adaptive(const RealType fMaxTimeStep)
{
//BField nie jest w tym prostym modelu uwiklane z P i moge potraktowac je jako funkcje t. Inaczej potrzebna bedzie wersja wektorowa.
	

	RealType fTimeStepSeconds =  MYears2Seconds(fMaxTimeStep);
	RealType fTimeStep        =  fMaxTimeStep;
	//RealType fAbsoluteEps     =  1e10*std::numeric_limits<RealType>::epsilon();
	RealType fAbsoluteEps     =  1e-3;//1e10*std::numeric_limits<RealType>::epsilon();

	RealType fA = m_PeriodDerivative(fPeriod, fBField);
	RealType fB = m_PeriodDerivative(fPeriod + fA*fTimeStepSeconds/2., m_NextBField(fBField, fTimeStep/2.));
	RealType fC = m_PeriodDerivative(fPeriod + fB*fTimeStepSeconds/2., m_NextBField(fBField, fTimeStep/2.));
	RealType fD = m_PeriodDerivative(fPeriod + fC*fTimeStepSeconds   , m_NextBField(fBField, fTimeStep));
	RealType fDPeriodNext = (fA + 2.*fB + 2.*fC + fD)*fTimeStepSeconds/6.;
	RealType fDPeriodPrevious(0.),fEps(0.),fError(0.);
	int i(1);

	do
	{
		fDPeriodPrevious = fDPeriodNext;
		fDPeriodNext     = 0.;

		i<<=1;//i*=2;
		fTimeStepSeconds /=2.;
		fTimeStep        /=2.;
		for(int j(0);j<i;j++)
		{
			fA = m_PeriodDerivative(fPeriod + fDPeriodNext                         , m_NextBField(fBField, fTimeStep*j)               );
			fB = m_PeriodDerivative(fPeriod + fDPeriodNext + fA*fTimeStepSeconds/2., m_NextBField(fBField, fTimeStep*j + fTimeStep/2.));
		 	fC = m_PeriodDerivative(fPeriod + fDPeriodNext + fB*fTimeStepSeconds/2., m_NextBField(fBField, fTimeStep*j + fTimeStep/2.));
			fD = m_PeriodDerivative(fPeriod + fDPeriodNext + fC*fTimeStepSeconds   , m_NextBField(fBField, fTimeStep*j + fTimeStep)   );
			fDPeriodNext += (fA + 2.*fB + 2.*fC + fD)*fTimeStepSeconds/6.;
		}
		fEps = fAbsoluteEps/fDPeriodNext;//Względny epsilon
		fError = pow(fTimeStep,5.)/fDPeriodNext;
		//cout<<i<<" "<<fDPeriodNext<<" "<<fDPeriodPrevious<<" "<<fEps<<" "<<fError<<endl;
	}
	while( !( ((fDPeriodPrevious-fEps) < fDPeriodNext) && ((fDPeriodPrevious+fEps) > fDPeriodNext) ) && !(fError < fEps));

	fPeriod += fDPeriodNext;
	//cout<<fBField;
	fBField = m_NextBField(fBField, fTimeStep*i);
	//cout<<" "<<fBField<<" "<<fTimeStep*i<<endl;
	fPeriodDot=m_PeriodDerivative(fPeriod,fBField);
}


RealType CPulsarPhysics::m_NextBField(const RealType fB, const RealType fTimeStep)
{
//cout<<m_pCfgPtr->Model.fBDecayDelta<<" "<<m_pCfgPtr->Model.fBDecayKappa<<" "<<fInitialBField<<" "<<fFinalBField<<" "<<fTimeStep+fAge<<endl;
//cout<<(fInitialBField - fFinalBField)*exp(-pow((fTimeStep+fAge)/m_pCfgPtr->Model.fBDecayDelta, m_pCfgPtr->Model.fBDecayKappa)) +fFinalBField<<endl;
  return (fB - fFinalBField)*exp(-pow((fTimeStep)/m_pCfgPtr->Model.fBDecayDelta, m_pCfgPtr->Model.fBDecayKappa)) +fFinalBField;
  //return (fB - fFinalBField)*exp(- fTimeStep/m_pCfgPtr->Model.fBDecayDelta ) +fFinalBField;
}

void CPulsarPhysics::m_EvolveBField(const RealType fTimeStep)
{
	fBField= m_NextBField(fBField, fTimeStep);
}






RealType CPulsarPhysics::m_fBFiledTimeDelatDistribution()
{
	return m_pRanPtr->Flat()*5. + 5.;
}

RealType CPulsarPhysics::m_fBFieldTimeDeltaGaussDistribution(const RealType fMean, const RealType fSigma)
{
	RealType fResult;
	do
	{
		fResult= m_pRanPtr->Gauss(fSigma, fMean);
	}
	while(fResult<=0.);
	return fResult;
}



CPulsarPhysics::CPulsarPhysics()
{
	m_pCfgPtr=CConfiguration::GetInstance();
	m_pRanPtr=m_pCfgPtr->pRan;


	
}


void CPulsarPhysics::Init()
{

	InitPulsarPhysics();


}


CPulsarPhysics::~CPulsarPhysics()
{

}


RealType CPulsarPhysics::m_FinalLogBFieldDistribution()
{
	
	RealType fParamLogBFieldMinMean(m_pCfgPtr->Model.fBFinalMean);
	RealType fParamLogBFieldMinSigma(m_pCfgPtr->Model.fBFinalSigma);

	return m_pRanPtr->Gauss(fParamLogBFieldMinSigma, fParamLogBFieldMinMean);
}

RealType CPulsarPhysics::m_InitLogBFieldDistribution()
{
	
	RealType fParamLogBFieldMean(12.65);//[G]
	RealType fParamLogBFieldSigma(0.55);//[G]
	
	if(m_pCfgPtr->MCMCParameterUse.bBInitMean) fParamLogBFieldMean=m_pCfgPtr->Model.fBInitMean;
	if(m_pCfgPtr->MCMCParameterUse.bBInitSigma) fParamLogBFieldSigma=m_pCfgPtr->Model.fBInitSigma;

	return m_pRanPtr->Gauss(fParamLogBFieldSigma, fParamLogBFieldMean);
	//return m_pRanPtr->Gauss(fParamLogBFieldSigma) + fParamLogBFieldMean;
}

RealType CPulsarPhysics::m_InitPeriodDistribution()
{
	RealType fParamPeriodMean(0.300);//[s]
	RealType fParamPeriodSigma(0.150);//[s]
	
	if(m_pCfgPtr->MCMCParameterUse.bPInitMean) fParamPeriodMean=m_pCfgPtr->Model.fPInitMean;
	if(m_pCfgPtr->MCMCParameterUse.bPInitSigma) fParamPeriodSigma=m_pCfgPtr->Model.fPInitSigma;

	//return m_pRanPtr->Gauss(fParamPeriodSigma, fParamPeriodMean);
	return m_pRanPtr->PositiveGauss(fParamPeriodSigma, fParamPeriodMean);
}


RealType CPulsarPhysics::m_InitLogPeriodDistribution()
{
	RealType fParamLogPeriodMean(-0.523);//odpowiada 0.3 dla lin
	RealType fParamLogPeriodSigma(0.2);
	
	if(m_pCfgPtr->MCMCParameterUse.bPInitMean) fParamLogPeriodMean=m_pCfgPtr->Model.fPInitMean;
	if(m_pCfgPtr->MCMCParameterUse.bPInitSigma) fParamLogPeriodSigma=m_pCfgPtr->Model.fPInitSigma;

	//return m_pRanPtr->Gauss(fParamPeriodSigma, fParamPeriodMean);
	return m_pRanPtr->Gauss(fParamLogPeriodSigma, fParamLogPeriodMean);
}

RealType CPulsarPhysics::m_PeriodDerivative(const RealType fP, const RealType fB)
{
	return fB*fB/(1.024e39*fP);
}

void CPulsarPhysics::InitPulsarPhysics()
{
	//fPeriod=pow(10.,m_InitLogPeriodDistribution());
	fPeriod=m_InitPeriodDistribution();
	fInitialBField=fBField=pow(10.,m_InitLogBFieldDistribution());
	fPeriodDot=m_PeriodDerivative(fPeriod,fBField);
	fFinalBField=pow(10.,m_FinalLogBFieldDistribution());
	fL400=0.;
	fS1400=0.;
	//m_fBFieldTimeDelta=m_fBFiledTimeDelatDistribution();

	fInitialPeriod=fPeriod;
	fInitialPeriodDot=fPeriodDot;	
	
	//m_fBFieldTimeDelta=m_fBFieldTimeDeltaGaussDistribution(m_pCfgPtr->Model.fBDecayTimeDeltaMean, m_pCfgPtr->Model.fBDecayTimeDeltaSigma);
}


#if 0
void CPulsarPhysics::m_ComputePPdotBField(const RealType fTimeStep)
{	
	//cout<<"Before: "<<fPeriod<<" "<<fBField<<" "<<fPeriodDot<<endl;
	const RealType fEta2(3.2e19 * 3.2e19);
	const RealType fKappa (1e6*1e7*M_PI);//seconds in 1Myr
	RealType fDB( fBField - m_fFinalBField );

#if DEBUG
	RealType NewfPeriod = sqrt( (fKappa/fEta2) * ( 2.*m_fFinalBField*m_fFinalBField*fTimeStep + 4.*fDB*m_fFinalBField*m_fBFieldTimeDelta*(1.-exp(-fTimeStep/m_fBFieldTimeDelta)) + fDB*fDB*m_fBFieldTimeDelta*(1.- exp( - 2.*fTimeStep/m_fBFieldTimeDelta)) )  + fPeriod*fPeriod);


	RealType NewfBField = m_NextBField(fBField, fTimeStep);
	RealType NewfPeriodDot = m_PeriodDerivative(fPeriod,fBField);
	if(NewfPeriod!=NewfPeriod)
	{
		cerr<<"Before:  "<<fTimeStep<<" "<<fPeriod<<" "<<fBField<<" "<<fPeriodDot<<endl;
		cerr<<"After:   "<<fTimeStep<<" "<<NewfPeriod<<" "<<NewfBField<<" "<<NewfPeriodDot<<endl;
	}

#endif
	fPeriod = sqrt( (fKappa/fEta2) * ( 2.*m_fFinalBField*m_fFinalBField*fTimeStep + 4.*fDB*m_fFinalBField*m_fBFieldTimeDelta*(1.-exp(-fTimeStep/m_fBFieldTimeDelta)) + fDB*fDB*m_fBFieldTimeDelta*(1.- exp( - 2.*fTimeStep/m_fBFieldTimeDelta)) )  + fPeriod*fPeriod);


	fBField = m_NextBField(fBField, fTimeStep);
	fPeriodDot = m_PeriodDerivative(fPeriod,fBField);

}
#endif


RealType CPulsarPhysics::m_RombergIntegralBFieldSquared(const RealType fTimeStep)
{
	const int nSize(4);
	RealType fTempBField(0.);
	double Tab[nSize];

	//Starting number of intervals
	int nNumber(4);
	double fH( (fTimeStep)/nNumber );


//Do poprawienia z BField
	//First integral
	//Tab[0] = (fH/2.) * ( m_NextBField(fBField,0.)*m_NextBField(fBField,0.) + m_NextBField(fBField, fTimeStep)*m_NextBField(fBField, fTimeStep));
	//And now the optimization. It's silly but the NextBField is almost 60% of the comp time (with dynamics switched off).
	fTempBField = m_NextBField(fBField,0.);
	Tab[0] += (fH/2.)*fTempBField*fTempBField;
	fTempBField = m_NextBField(fBField, fTimeStep);
	Tab[0] += (fH/2.)*fTempBField*fTempBField;


	for(int i(1); i<nNumber; i++)
	{
		//Tab[0] += fH*m_NextBField(fBField,0. + i*fH)*m_NextBField(fBField,0. + i*fH);
		fTempBField = m_NextBField(fBField,0. + i*fH);
		Tab[0] += fH*fTempBField*fTempBField;
	}

	//Rest of integrals
	for(int i(1);i<nSize;i++)
	{	
		Tab[i] = Tab[i-1]/2.;
		for(int j(0); j<nNumber; j++)
		{
			//Tab[i] += (fH/2.) * m_NextBField(fBField,fH/2. + j*fH)*m_NextBField(fBField,fH/2. + j*fH);
			fTempBField = m_NextBField(fBField,fH/2. + j*fH);
			Tab[i] += (fH/2.) * fTempBField*fTempBField; 
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



RealType CPulsarPhysics::m_PeriodFunction(const RealType fTimeStep)
{
	const RealType fEta2(3.2e19 * 3.2e19);//const terms from the lighthouse model squared
	const RealType fKappa (1e6*1e7*M_PI);//seconds in 1Myr

	return 	sqrt((fKappa/fEta2) * \
			( fFinalBField*fFinalBField*fTimeStep + \
			(1.-exp(-2.*fTimeStep/m_pCfgPtr->Model.fBDecayDelta))*(fBField*fBField - 2.*fBField*fFinalBField + fFinalBField*fFinalBField)*(m_pCfgPtr->Model.fBDecayDelta/2.) + \
			(1.-exp(-fTimeStep/m_pCfgPtr->Model.fBDecayDelta))*(fFinalBField*fBField-fFinalBField*fFinalBField)*2.*m_pCfgPtr->Model.fBDecayDelta ) \
		+fInitialPeriod*fInitialPeriod );
}




void CPulsarPhysics::m_ComputePPdotBField(const RealType fTimeStep)
{	
	const RealType fEta2(3.2e19 * 3.2e19);
	const RealType fKappa (1e6*1e7*M_PI);//seconds in 1Myr

	//fPeriod = sqrt((fKappa/fEta2)*m_RombergIntegralBFieldSquared(fTimeStep) + fPeriod*fPeriod);

	fPeriod = m_PeriodFunction(fTimeStep);

	//cout<<"ra: "<<fPeriod<<" "<<fTemp<<" "<<fPeriod-fTemp<<" "<<fTimeStep<<" "<<fInitialPeriod<<" "<<fInitialBField<<" "<<fFinalBField<< endl;


	fBField = m_NextBField(fBField, fTimeStep);
	fPeriodDot = m_PeriodDerivative(fPeriod,fBField);

}


void CPulsarPhysicTester::InitPeriod(const int nNum, const char* FileName)
{
	RealType LogP;
	ofstream out(FileName);
	for(int i(nNum);i-- >0;)
	{
		InitPulsarPhysics();
		LogP=log10(this->fPeriod);
		//LogP=m_InitLogPeriodDistribution();
		out<<LogP<<" "<<pow(10.,LogP)<<endl;
	}
	out.close();
}


void CPulsarPhysicTester::FinalBField(const int nNum, const char* FileName)
{
	RealType LogB;
	ofstream out(FileName);
	for(int i(nNum);i-- >0;)
	{
		InitPulsarPhysics();
		LogB=log10(this->fFinalBField);
		//LogB=m_FinalLogBFieldDistribution();
		out<<LogB<<" "<<pow(10.,LogB)<<endl;
	}
	out.close();
}


void CPulsarPhysicTester::InitBField(const int nNum, const char* FileName)
{
	RealType LogB;
	ofstream out(FileName);
	for(int i(nNum);i-- >0;)
	{
		InitPulsarPhysics();
		LogB=log10(this->fBField);
		//LogB=m_InitLogBFieldDistribution();
		out<<LogB<<" "<<pow(10.,LogB)<<endl;
	}
	out.close();
}
