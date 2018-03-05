#include<RadioDetector.hpp>

#include<cmath>
#include<Debug.hpp>
#include<fstream>
#include<string>

using namespace std;

#include<iostream>


RealType CRadioDetector::m_Flux(const CPulsar *pPulsar, const CSurvey* pSurvey) const
{
	return m_Flux400toFreq(pPulsar, pSurvey->fFrequency);
}

RealType CRadioDetector::m_Flux400toFreq(const CPulsar *pPulsar, const RealType fFrequency) const
{
	//RealType S400 = pPulsar->fL400 / (pPulsar->GalacticPos.fDist*pPulsar->GalacticPos.fDist);
	//praca Osłowski, Bulik ->
	//const RealType fKsi(-1.8);
	//return S400*pow( fFrequency/400., fKsi );

	return pPulsar->RadioFlux(fFrequency);
}



bool CRadioDetector::CheckParkes(const CPulsar *pPulsar) const
{

	return IsVisible(pPulsar, ParkesMultiBeam);

}


bool CRadioDetector::IsVisible(const CPulsar *pPulsar, const ESurveyID ESurvey) const
{
    if(static_cast<int>(ESurvey)>gSurveyNum) ERROR("Something is off...");
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
    CSurvey* pSurvey=&(m_pSurvey[static_cast<int>(ESurvey)]);

    //cout<<ESurvey<<" "<<pSurvey->Name<<endl;

	if(m_Flux(pPulsar,pSurvey)>=m_SMinimal(pPulsar, pSurvey))
	{
		if(pCfgPtr->Sim.bSurveyGeometry ) return pSurvey->IsCovered(pPulsar->GalacticPosDeg.fL,pPulsar->GalacticPosDeg.fB);
		else return true;
	}

	return false;
}




CRadioDetector::CRadioDetector()
{
	m_SNRTreshold=10.;
	m_pTemperatureSky = new short int[360*180];
	m_ReadTemperatureSkyData("input/tsky1.ascii");
	m_nSurveyNum = gSurveyNum;//where is it located?
	m_pSurvey = new CSurvey[m_nSurveyNum];
	for(int i(0);i<gSurveyNum;i++)
	{
		m_pSurvey[i].Init(static_cast<ESurveyID>(i));
	}

}


CRadioDetector::~CRadioDetector()
{
	delete [] m_pTemperatureSky;
	delete [] m_pSurvey;
}



RealType CRadioDetector::m_TauScatter(const CPulsar *pPulsar, const CSurvey* pSurvey) const
{
	RealType fTauScatter;
   	//Ramesh Bhat, Cordes, Lazio... 2004 equation 7:
   	fTauScatter = -6.46 + 0.154 * log10(pPulsar->fDM) + 1.07 * log10(pPulsar->fDM) * log10(pPulsar->fDM) - 3.86 * log10(pSurvey->fFrequency * 1e-3);
   	fTauScatter = pow(10., fTauScatter); // ms
   	fTauScatter = fTauScatter * 1e-3; // s
   	return fTauScatter;
}

RealType CRadioDetector::m_TemperatureSky(const Galactic PosDeg, const CSurvey* pSurvey) const
{
	RealType fWavelength = gCLight/(pSurvey->fFrequency*1e6);
	RealType fTSky400    = m_TemperatureSkyTab(PosDeg.fL,PosDeg.fB);
	RealType fLambda408  = gCLight/(408.*1e6);/* wavelength corresponding to 408MHz, m */
	return  fTSky400 * pow(fWavelength/fLambda408, 2.8);
}


// returns We, effective width of the pulsar (observed in a survey)
RealType CRadioDetector::m_WEffective(const CPulsar *pPulsar, const CSurvey* pSurvey) const
{
	RealType fTauScatter = m_TauScatter(pPulsar,pSurvey);
	return sqrt(pow(m_WIntrinsic(pPulsar),2) + (1. + pow(pPulsar->fDM/pSurvey->fDM0,2)) * (pSurvey->fSamplingTime * pSurvey->fSamplingTime) + fTauScatter * fTauScatter);
}

RealType CRadioDetector::m_WIntrinsic(const CPulsar *pPulsar) const
{  
   	//return pPulsar->fFWHM * pPulsar->fPeriod; // Width in [s]
	return 0.05 * pPulsar->fPeriod;
}

RealType CRadioDetector::SMinimal(const CPulsar *pPulsar, const ESurveyID ESurvey) const
{
    if(static_cast<int>(ESurvey)>gSurveyNum) ERROR("Something is off...");
	return	m_SMinimal(pPulsar, &(m_pSurvey[static_cast<int>(ESurvey)]) );
}


RealType CRadioDetector::m_SMinimal(const CPulsar *pPulsar, const CSurvey* pSurvey) const
{

	RealType fSMinimal;
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	RealType fDetectionTreshold(pCfgPtr->Sim.fDetectionTreshold);	
//	const RealType cfDetectionTreshold(10.);   // Detection threshold.
	int nNPol(2); //Liczba obserwowanych polaryzacji

   	// Work out effective width of the pulsar.
	RealType fWEffective =  m_WEffective(pPulsar,pSurvey);
	if(fWEffective!=fWEffective)
	{
		WARNING("Problem computing effective width");
		CPulsar Dupa=*pPulsar;
		Dupa.Project();
		fWEffective =  m_WEffective(pPulsar,pSurvey);
	}

   	if (fWEffective >= pPulsar->fPeriod)
	{
      /* The telescope cannot 'see' the pulsar. */
      /* Pulsar should never be detected. */
      		fSMinimal = 9.99e15;
   	}
   	else
   	{  // vals are surv dep
      		fSMinimal = sqrt(M_PI/2.) * sqrt(fWEffective/(pPulsar->fPeriod-fWEffective))*fDetectionTreshold*(pSurvey->fTemperatureReceiver + m_TemperatureSky(pPulsar->GalacticPosDeg, pSurvey) ) / ( pSurvey->fGain * sqrt(nNPol*pSurvey->fReceiverBandwidth  * pSurvey->fIntegrationTime ));
   	}
	if(fSMinimal!=fSMinimal)
	{
		cout<<"DDD IsLost: "<<pPulsar->bIsLost<<endl;
		cout<<"DDD SMin: "<<fSMinimal<<"[Jy]"<<endl;
		cout<<"DDD P: "<<pPulsar->fPeriod<<" PDot: "<<pPulsar->fPeriodDot<<" DM: "<<pPulsar->fDM <<" Weff: "<<fWEffective  <<endl;	
		cout<<"DDD X: "<<pPulsar->Position.fX<<" Y: "<<pPulsar->Position.fY<<" Z: "<<pPulsar->Position.fZ<<endl;	

	}
   	return fSMinimal;
   
}


// returns S/N for pulsar with flux S

RealType CRadioDetector::m_SNR(const CPulsar *pPulsar, const CSurvey* pSurvey) const
{
	RealType fWEffective =  m_WEffective(pPulsar,pSurvey);
	RealType fSNR(0.);
	RealType fFlux=m_Flux(pPulsar,pSurvey);
  	int nNPol(2);

	if (fWEffective < pPulsar->fPeriod)
  	{
	    fSNR = fFlux / sqrt(M_PI/2.) / sqrt(fWEffective/(pPulsar->fPeriod-fWEffective)) / 
			(pSurvey->fTemperatureReceiver + m_TemperatureSky(pPulsar->GalacticPosDeg,pSurvey) ) * (pSurvey->fGain * sqrt(nNPol*pSurvey->fReceiverBandwidth * pSurvey->fIntegrationTime ));
	}

  	return fSNR;
}



void CRadioDetector::m_ReadTemperatureSkyData(const char * Filename)
{
	std::ifstream in(Filename);
	if(!in.good())
	{
		string msg("Could not open SkyTemprature datafile ");
		msg+=Filename;
		ERROR(msg);
	}

	for(int i(0);i<360*180;i++)
	{
		in>>m_pTemperatureSky[i];
	}
	in.close();

}

RealType CRadioDetector::m_TemperatureSkyTab(const RealType fLDeg, const RealType fBDeg) const
{
	if( fLDeg < 0. || fLDeg >= 360. || fBDeg < -90. || fBDeg > 90.)
	{
		cerr<<fLDeg<<" "<<fBDeg<<endl;
		ERROR("Galactic coordinates out of range");
	}
	//2Deg
	//RealType fLDeg(fL*180./M_PI);
	//RealType fBDeg(fB*180./M_PI);
	//mapa jest z rozdzielczosci co 1deg, wyposrodkowane w calkowitych
	//za sky1.f:
	int nX=static_cast<int>(fLDeg+0.5);
	int nY=static_cast<int>(fBDeg+90.5);
	if(nX>=360) nX=359;//wydaje mi sie, ze raczej powinno byc 0 //FIXME sprawdzic w papierze
	if(nX>=180) nY=179;//wydaje mi sie, ze raczej powinno byc 0 //FIXME sprawdzic w papierze
	return 0.1*m_pTemperatureSky[nX + 360*nY];
}

void CRadioDetector::Print4StefansTest(const CPulsar *pPulsar) const
{
	if(! pPulsar->bIsLost)
	{
CSurvey Parkes(Parkes70);


	Cartesian PosRot = RotateVectorXYPlane(pPulsar->Position, M_PI/2.);

cout<<"PosXYZ: "<<PosRot.fX<<" "<<PosRot.fY<<" "<<PosRot.fZ<<" ";
cout<<"VelXYZ: "<<pPulsar->Velocity.fX<<" "<<pPulsar->Velocity.fY<<" "<<pPulsar->Velocity.fZ<<" ";
cout<<"P: "<<pPulsar->fPeriod<<" Pdot: "<<pPulsar->fPeriodDot<<" "<<0.05<<" "<<0<<endl;
cout<<"b: "<<pPulsar->GalacticPosDeg.fB<<" l: "<<pPulsar->GalacticPosDeg.fL<<" d: "<<pPulsar->GalacticPosDeg.fDist<<" DM: ";
cout<<pPulsar->fDM<<" DM0: "<<Parkes.fDM0<<" Tsky: "<<m_TemperatureSky(pPulsar->GalacticPosDeg, &Parkes)<<" Tau "<<m_TauScatter(pPulsar,&Parkes)<<" We "<<m_WEffective(pPulsar,&Parkes)<<" Smin "<<m_SMinimal(pPulsar, &Parkes)<<" S400: "<<m_Flux(pPulsar,&Parkes)<<" Vis:"<<CheckParkes(pPulsar)<<endl;






	}
} 

