#include<Pulsar.hpp>
#include<iostream>
#include<cmath>
#include<limits>
#include<sstream>
using namespace std;


//should be for the L400
RealType CPulsar::LEfficiency() const
{
	if(m_pCfgPtr->Sim.bLEfficiency)
	{
		return 1.;
	}

	return 1.;
}

RealType CPulsar::LogLCorr() const
{
	const RealType fSigmaCorr(0.8);
	if(m_pCfgPtr->Sim.bLogLCorr)
	{
		return m_pRanPtr->Gauss(fSigmaCorr);
	}

	return 0.;
}

CPulsar::CPulsar()
{
	m_pCfgPtr=CConfiguration::GetInstance();
	m_pRanPtr=m_pCfgPtr->pRan;
	m_pBS = static_cast<SBasicStar*>(this);
	//m_pBS = dynamic_cast<SBasicStar*>(this);
}

CPulsar::~CPulsar()
{
}

CPulsar::CPulsar(const SBasicStar& Star)
{
	CPulsar();
	*m_pBS=Star;
}

CPulsar& CPulsar::operator=(const SBasicStar& BS)
{
	*m_pBS=BS;
	return *this;
}

#if 0
//TODO FIXME 
//zmienic na operator <<
//i przeniesc do BasiStar
void CPulsar::Print() const
{
	PrintStarInfo(std::cout);
}


#endif


bool CPulsar::m_CheckDistance()
{
	RealType fDistanceFromGalCenter(sqrt(Position.fX*Position.fX + Position.fY*Position.fY + Position.fZ*Position.fZ));
	if(fDistanceFromGalCenter > m_pCfgPtr->Sim.fMaxDist)
	{
		return false;
	}

	return true;	

}


void CPulsar::Project()
{
	
	if(m_pCfgPtr->Sim.bDynamics)
	{
		GalacticPos=m_Project();
		GalacticPosDeg=GalacticPos;
		GalacticPosDeg.fL *= (180./M_PI);
		GalacticPosDeg.fB *= (180./M_PI);
		if( m_CheckDistance() && m_Check4Relevance() )
		{
			fDM=DispersionMeasure();
			bInsideGalaxy = true;
		}
		else
		{
			bInsideGalaxy = false;
		}
	}



	if(m_pCfgPtr->Sim.bPhysics) 
	{
		RealType fDeath(DeathLinesProbability());
		bool bLife = m_pRanPtr->Flat() < fDeath;
		fProbability *= fDeath;
		bBeaming = false;
		if(bLife)
		{
			bBeaming = m_CheckBeaming();
		//bIsLost = !(  bInsideGalaxy && bBeaming && CheckProbability() );
		}	
		bIsLost = !(  bInsideGalaxy && bBeaming && bLife );
	}
	else
	{
		bIsLost = !bInsideGalaxy;
	}

	if(m_pCfgPtr->Sim.bObservations && !bIsLost)
	{
		fL400= pow(10.,LogRadioLuminosity());///< [mJy*kpc^2]
//		fL400= pow(10.,LogL400());///< [mJy*kpc^2]
		fS400 = RadioFlux400();
		fS1400 = RadioFlux(1400.); ///< [mJy]
	}
	//cout<<"Project "<<fS1400<<" "<<fProbability<<" "<<bBeaming<<" "<<bIsLost<<endl;
}

void CPulsar::EvolveDynamics(const RealType fMaxTimeStep)
{
	if(m_pCfgPtr->Sim.bDynamics) CPulsarDynamics::Evolve(fMaxTimeStep);
}
void CPulsar::EvolvePhysics(const RealType fMaxTimeStep)
{
	if(m_pCfgPtr->Sim.bPhysics)  CPulsarPhysics::Evolve(fMaxTimeStep);
}

void CPulsar::Evolve(const RealType fMaxTimeStep)
{

	
	if(m_pCfgPtr->Sim.bDynamics) CPulsarDynamics::Evolve(fMaxTimeStep);
	if(m_pCfgPtr->Sim.bPhysics)  CPulsarPhysics::Evolve(fMaxTimeStep);
}


void CPulsar::Init(const RealType fPositionStartTime)
{
	fSpectralIndex=SpectralIndexDistribution();
	fProbability=1.;
	fMass=1.4;

	if(m_pCfgPtr->Sim.bDynamics) CPulsarDynamics::Init(fPositionStartTime);
	if(m_pCfgPtr->Sim.bPhysics)  CPulsarPhysics::Init();

	if(!m_pCfgPtr->Sim.bDynamicsComputed)
	{
		bBeaming = true;
		bInsideGalaxy = true;	

	}
}




RealType CPulsar::SpectralIndexDistribution() const
{
	//Maron at al. 2000
	//return -1.8;
	//Bates at al 2013
	return -1.4;
}


bool CPulsar::CheckProbability() const
{
	const RealType fAbsoluteEps = 1e6*std::numeric_limits<RealType>::epsilon();
	return (fProbability>fAbsoluteEps);
}

RealType CPulsar::DeathLinesProbability() const
{
	//return 1.*AboveDeathLines();
	return LifeProbability(log10(fPeriod), log10(fPeriodDot), m_pCfgPtr->Model.fDeathPsi);
}

//starsz funkcja
bool CPulsar::AboveDeathLines() const
{
	return AboveDeathLines(fPeriod, fPeriodDot);
}

bool CPulsar::AboveDeathLines(const RealType fP, const RealType fPDot) const
{
//	const RealType fMi = - 0.2;
	const RealType fMi = 0.;
	bool A = log10(fPDot) >= 3.29 *log10(fP)-16.55;
	bool B = log10(fPDot) >= 0.92*log10(fP)-18.65 + fMi;
//cout<<"AboveDeathLines "<< log10(fPeriodDot)<<" "<<3.29 *log10(fPeriod)-16.55<<endl;
//cout<<"AboveDeathLines "<< log10(fPeriodDot)<<" "<<0.92 *log10(fPeriod)-18.65+fMi<<endl;
//cout<<"AboveDeathLines "<<A<<" "<<B<<" "<<(A && B)<<endl;
	return (A && B);
}


RealType CPulsar::LogRadioLuminosity()
{
	//FIXME wyprowadzić 
	//if(m_pCfgPtr->Model.bLPowerLaw) return LogL400PowerLaw();
	//if(m_pCfgPtr->Model.bLRandom) return LogL1400Random()*pow(400./1400.,fSpectralIndex);;
	//WARNING("No Luminosity-model flag specified. Narayan, Ostriker model used.");
	//
	//
	RealType fLogL400(0.);
	switch(m_pCfgPtr->Sim.LuminosityModel)
	{	
	case Narayan:
		fLogL400 = LogL400();
		break;
	case Proszynski:
		fLogL400 = log10( L400Proszynski() );
		break;
	case PowerLaw:
		fLogL400 = LogL400PowerLaw();
		break;
	case Erot:
		fLogL400 = LogL400Erot();
		break;
	case Random:
		fLogL400 = log10(RadioLuminosity400(pow(10.,LogL1400Random()), 1400.) );
		break;
	default:
		ERROR("No luminosity model specified.");
		break;
	}

	return LEfficiency()*( fLogL400 + LogLCorr() );
}


RealType CPulsar::LogL1400Random() const
{
	const RealType fSigma(0.9);
	const RealType fMi(-1.1);

	return m_pRanPtr->Gauss(fSigma,fMi);
}


RealType CPulsar::LogL400PowerLaw() const
{
	return log10( m_pCfgPtr->Model.fLPowerLawGamma * pow(fPeriod,m_pCfgPtr->Model.fLPowerLawAlpha) * pow(fPeriodDot*1e15, m_pCfgPtr->Model.fLPowerLawBeta));


}


//Narayan, Ostriker 1990
RealType CPulsar::LogL400()
{
	RealType fResult =(1./3.) * ( log10(fPeriodDot*1e15) - log10(fPeriod*fPeriod*fPeriod) ) + 1.635; //log([mJy*kpc^2])
	if(fResult!=fResult)
	{
		cerr<<fPeriod<<" "<<fPeriodDot<<" "<<log10(fPeriodDot)<<" "<<log10(fPeriod)<<endl;
	}

	return fResult;
}



#if 0
RealType LogL1400()
{
	RealType fL0;
	
	RealType fSigmaLCorr;
	RealType fLCorr = m_pRanPtr->Gauss(fSigmaLCorr);
	return log10(fL0*pow(fPeriod,fEpsilonP)*pow(fPeriodDot/(1e15),fEpsilonP ) ) + fLCorr;
}
#endif

RealType CPulsar::DispersionMeasure()
{
	float dummy(0.);
	RealType fDM;

	if(m_pCfgPtr->Sim.bComputeDMDiff)
	{
		return m_pCfgPtr->pDM->GetDM(GalacticPosDeg.fL, GalacticPosDeg.fB, GalacticPosDeg.fDist) - fort_dmdsm(GalacticPos.fL, GalacticPos.fB, -1, dummy, GalacticPos.fDist);
	}


	if(m_pCfgPtr->Sim.bUsePrecomputedDMData)
	{
		fDM=m_pCfgPtr->pDM->GetDM(GalacticPosDeg.fL, GalacticPosDeg.fB, GalacticPosDeg.fDist);
	}
	else 
	{
		fDM=fort_dmdsm(GalacticPos.fL, GalacticPos.fB, -1, dummy, GalacticPos.fDist);
	}

	return fDM;
}


RealType CPulsar::L400Proszynski() const
{
	return pow(10.,7.01)*pow(fPeriod,-1.04)*pow(fPeriodDot,0.35);
}
RealType CPulsar::RadioLuminosity400(const RealType fLuminosity,const RealType fFrequency ) const
{
	return fLuminosity*pow( fFrequency/400., fSpectralIndex );	
}

RealType CPulsar::RadioFlux400() const
{
	return fL400 / (GalacticPos.fDist*GalacticPos.fDist);
}

RealType CPulsar::RadioFlux(const RealType fFrequency) const
{
//	RealType S400 = fL400 / (GalacticPos.fDist*GalacticPos.fDist);
	return fS400*pow( fFrequency/400., fSpectralIndex );		
}
 ///result in [Jy], freq in [MHz]
	
bool CPulsar::m_CheckBeaming()
{
	if( m_pRanPtr->Flat() <= BeamingFactor(fPeriod) ) return true;
	return false;
}
RealType CPulsar::BeamingFactor() const
{
	return BeamingFactor(fPeriod);
}
RealType CPulsar::BeamingFactor(const RealType fP) const
{

	return (9.*pow(log10(fP/10.),2.)+3.)/100.;
} 


bool CPulsar::m_Check4Relevance() const
{
	if(m_pCfgPtr->Sim.bCheck4Relevance)
	{
		if(!m_pCfgPtr->SurveysCuts.bInitiated) ERROR("bCheck4Relevance is true but SurveysCuts are not initiated");

		for(int i(0); i<m_pCfgPtr->SurveySelection.nNum; ++i)
		{
			bool bDist = m_pCfgPtr->SurveysCuts.pMaxSurveyDist[i] >GalacticPos.fDist;
			bool bCoordinates = CSurvey::IsCoveredInSurvey(GalacticPosDeg.fL, GalacticPosDeg.fB, m_pCfgPtr->SurveySelection.pID[i]);
			if(bDist && bCoordinates) return true;
			else return false;
		}
	}
	return true;
}


RealType CPulsar::LogL400Erot() const
{
	return log10( m_pCfgPtr->Model.fLPowerLawGamma * pow( (pow(fPeriodDot*1e15,1./3.) / fPeriod), m_pCfgPtr->Model.fLPowerLawAlpha) );
}




RealType  CPulsar::LogPFromDeathLine(const RealType fLogPDot) const
{
	if (fLogPDot > -19.7424)
	{
		return (fLogPDot+16.55)/3.29;
	}
	else
	{
		//const RealType fMi(-0.2);
		const RealType fMi(0.);
		return (fLogPDot+18.65-fMi)/0.92;
	}
}

RealType  CPulsar::LifeProbability(const RealType fLogP, const RealType fLogPDot, const RealType fPhi) const
{
	if(m_pCfgPtr->Sim.bNoDeathLine)
	{
		return 1.;
	}
	return -atan( (fLogP-LogPFromDeathLine(fLogPDot) )/fPhi)/M_PI+0.5;
}
   

void  CPulsarTester::DeathTester() const
{    
	ofstream out("DeathArea.mtx");
	const int nRes(100);
	for(RealType fLogP(-2.); fLogP<=1.; fLogP+=3./nRes)
	{
		for(RealType fLogPDot(-25.); fLogPDot<=-10.; fLogPDot+=15./nRes)
		{
			out<<LifeProbability(fLogP, fLogPDot, m_pCfgPtr->Model.fDeathPsi)<<" ";
		}
		out<<endl;
	}
	out.close();
}



void  CPulsarTester::LErotTester(const RealType fGamma, const RealType fAlpha)
{

	cout<<fGamma<<" "<<fAlpha<<endl;
	m_pCfgPtr->Model.fLPowerLawGamma=fGamma;
	m_pCfgPtr->Model.fLPowerLawAlpha=fAlpha;
	ofstream out("LErot.mtx");
	const int nRes(100);
	for(RealType fLogP(-2.); fLogP<=1.; fLogP+=3./nRes)
	{
		for(RealType fLogPDot(-25.); fLogPDot<=-10.; fLogPDot+=15./nRes)
		{
			fPeriod = pow(10.,fLogP);
			fPeriodDot = pow(10., fLogPDot);
			
			out<<LogL400Erot()<<" ";
		}
		out<<endl;
	}
	out.close();



}

std::string CPulsar::AsCatalogueEntry(const int nRefNumber) const
{
	std::stringstream Out;
	std::string Nothing("nan");

	Out << nRefNumber <<" ";//1 numer porzadkowy
	Out << Nothing <<" ";//2 nazwa
	Out << Nothing <<" ";//3 ref
	Out << Nothing <<" ";//4 jeszcze raz nazwa
	Out << Nothing <<" ";//5 ref

	Out << GalacticPosDeg.fL <<" ";//6 GalL
	Out << GalacticPosDeg.fB <<" ";//6 GalB
	Out << fPeriod<<" ";//8 P0

	Out << Nothing <<" ";//9 dP0
	Out << Nothing <<" ";//10 ref

	Out << fPeriodDot<<" ";//8 P0


	Out << Nothing <<" ";//12 dP1
	Out << Nothing <<" ";//13 ref

	

	Out << fDM<<" ";//14 DM
	Out << Nothing <<" ";//15 dDM
	Out << Nothing <<" ";//16 ref
	Out << Nothing <<" ";//17 RM
	Out << Nothing <<" ";//18 dRM
	Out << Nothing <<" ";//19 ref
	Out << Nothing <<" ";//20 binary
	Out << Nothing <<" ";//21 ref

	Out << GalacticPosDeg.fDist<<" ";//22 Dist
	Out << GalacticPosDeg.fDist<<" ";//23 DistDM

	Out << Nothing <<" ";//24 ref
	
	Out <<"pksmb ";//25 Survey
	Out <<"40000 ";//26 surevy mask
	
	Out << Nothing <<" ";//27 Radio Flux @ 400Mhz [mJy]
	Out << Nothing <<" ";//28 dS400
	Out << Nothing <<" ";//29 ref


	Out << fS1400<<" ";//30 Radio Flux @ 400Mhz [mJy]
	Out << Nothing <<" ";//31 dS1400
	Out << Nothing <<" ";//32 ref

	Out << fAge<<" ";//33 Age [yr]
	Out << fBField<<" ";//34 BSurf[G]
	Out << Nothing <<" ";//35 P1_I
	Out << Nothing <<" ";//36 Age_I
	Out << Nothing <<" ";//37 Bsurf_I
	Out << endl;

	return Out.str();
}

RealType CPulsar::ETot() 
{
	
	RealType fVTot( Length(Velocity) );
	RealType fEKin( fMass*fVTot*fVTot*0.5 );
	RealType fEPot( fMass*GalacticPotential(Position) );

	return fEKin - fEPot;
}

void CPulsarTester::CatModelLTest() 
{
    RealType fCatDist, fCatS1400;
    for(int i(0);i<m_pCfgPtr->pATNF->nTotalNum;i++)
    {
        if(m_pCfgPtr->pATNF->pbUsed[i])
        {


            fCatDist = m_pCfgPtr->pATNF->pPulsar[i].fDist;
            fCatS1400 = m_pCfgPtr->pATNF->pPulsar[i].fRadioFlux1400;

            fPeriod=m_pCfgPtr->pATNF->pPulsar[i].fP;
            fPeriodDot=m_pCfgPtr->pATNF->pPulsar[i].fPdot;
            GalacticPos.fDist = fCatDist;
            fL400 = pow(10.,LogRadioLuminosity());
            fS400 = RadioFlux400();
            fS1400 = RadioFlux(1400.);

            cout<<fCatS1400<<" "<<fS1400<<endl;
        }
    }

}
