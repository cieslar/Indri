#include<Survey.hpp>

#include<Definitions.hpp>
#include<Debug.hpp>
#include<cmath>
#include<Coordinates.hpp>
#include<cstring>
#include<iostream>
#include<strlcpy.hpp>

using namespace std;

CSurvey::CSurvey(ESurveyID SrvID) { Init(SrvID); }

void CSurvey::Init(const ESurveyID SrvID) {
    ID = SrvID;
    if (SrvID==Parkes70) {
        strlcpy(Name, "Parkes 70 cm", nNameLength);
        fTemperatureReceiver = 60.;
        fReceiverBandwidth = 32.;
        nChannelsNum = 256;
        fFrequency = 436.;
        fSamplingTime = 0.0003;
        fGain = 0.65;
        fIntegrationTime = 157;
        IsCovered = &SkyCoverageParkes70;
        fDM0 = m_DM0();
        return;
    }

    if (SrvID==ParkesMultiBeam) {
        strlcpy(Name, "Parkes Multi Beam", nNameLength);
        fTemperatureReceiver = 24.;
        fReceiverBandwidth = 288.;
        nChannelsNum = 96;
        // fFrequency=1374.;
        fFrequency = 1400.; // Freq the same as in ATNF
        fSamplingTime = 0.00025;
        fGain = 0.65;
        fIntegrationTime = 2100.;
        IsCovered = &SkyCoverageParkesMB;
        fDM0 = m_DM0();
        return;
    }

    if (SrvID==SwinIntermLat) {
        strlcpy(Name, "Swin Interm Lat", nNameLength);
        fTemperatureReceiver = 24.;
        fReceiverBandwidth = 288.;
        nChannelsNum = 96;
        fFrequency = 1374.;
        fSamplingTime = 0.000125;
        fGain = 0.65;
        fIntegrationTime = 265.;
        IsCovered = &SkyCoverageSwinIL;
        fDM0 = m_DM0();
        return;
    }

    if (SrvID==SwinExtended) {
        strlcpy(Name, "Swin Extended", nNameLength);
        fTemperatureReceiver = 24.;
        fReceiverBandwidth = 288.;
        nChannelsNum = 96;
        fFrequency = 1374.;
        fSamplingTime = 0.000125;
        fGain = 0.65;
        fIntegrationTime = 265.;
        IsCovered = &SkyCoverageSwinExt;
        fDM0 = m_DM0();
        return;
    }

    if (SrvID==Burgay) {
        strlcpy(Name, "Burgay et al", nNameLength);
        fTemperatureReceiver = 24.;
        fReceiverBandwidth = 288.;
        nChannelsNum = 96;
        fFrequency = 1374.;
        fSamplingTime = 0.000125;
        fGain = 0.65;
        fIntegrationTime = 265.;
        IsCovered = &SkyCoverageBurgay;
        fDM0 = m_DM0();
        return;
    }

    if (SrvID==GMRT) {
        strlcpy(Name, "GMRT", nNameLength);
        fTemperatureReceiver = 90.;
        fReceiverBandwidth = 32.;
        nChannelsNum = 512;
        fFrequency = 610.;
        fSamplingTime = 0.000125;
        fGain = 0.65 * 15. * (45.0 / 64.0 * 45.0 / 64.0);
        fIntegrationTime = 300.;
        IsCovered = &SkyCoverageGMRT;
        fDM0 = m_DM0();
        return;
    }

    if (SrvID==ParkesMultiBeamALLSKY) {
        strlcpy(Name, "Parkes Multi Beam ALLSKY", nNameLength);
        fTemperatureReceiver = 24.;
        fReceiverBandwidth = 288.;
        nChannelsNum = 96;
        fFrequency = 1374.;
        fSamplingTime = 0.00025;
        fGain = 0.65;
        fIntegrationTime = 2100.;
        IsCovered = &SkyCoverageAllSky;
        fDM0 = m_DM0();
        return;
    }

    if (SrvID==ParkesMultiBeamPart) {
        strlcpy(Name, "Parkes Multi Beam part", nNameLength);
        fTemperatureReceiver = 24.;
        fReceiverBandwidth = 288.;
        nChannelsNum = 96;
        fFrequency = 1374.;
        fSamplingTime = 0.00025;
        fGain = 0.65;
        fIntegrationTime = 2100.;
        IsCovered = &SkyCoveragePartSky;
        fDM0 = m_DM0();
        return;
    }

    if (SrvID == SKA1Mid) {
        strlcpy(Name, "SKA-1 Mid", nNameLength);
        fTemperatureReceiver = 30.;
        fReceiverBandwidth = 300.;
        nChannelsNum = 4096;
        // fFrequency=1374.;
        fFrequency = 1400.; // Freq the same as in ATNF
        fSamplingTime = 0.000064;
        fGain = 15.;
        fIntegrationTime = 900.; //[s]
        IsCovered = &SkyCoverageSKA1Mid;
        fDM0 = m_DM0();
        return;
        // cout<<"SKA-1 Mid DDM: "<<fDM0<<endl;
        // PrintStacktrace();
    }

    if (SrvID==SKA1MidLC) {
        strlcpy(Name, "SKA-1 Mid LowCost", nNameLength);
        fTemperatureReceiver = 30.;
        fReceiverBandwidth = 300.;
        nChannelsNum = 1024;
        // fFrequency=1374.;
        fFrequency = 1400.; // Freq the same as in ATNF
        fSamplingTime = 0.000064;
        fGain = 2.;
        fIntegrationTime = 900.; //[s]
        IsCovered = &SkyCoverageSKA1Mid;
        fDM0 = m_DM0();
        return;
        // cout<<"SKA-1 Mid DDM: "<<fDM0<<endl;
        // PrintStacktrace();
    }

    ERROR("No such survey");
    return;
}

bool CSurvey::IsCoveredInSurvey(const RealType fLDeg, const RealType fBDeg,
                                const ESurveyID SrvID) {
    if (SrvID==ParkesMultiBeam)
        return SkyCoverageParkesMB(fLDeg, fBDeg);
    if (SrvID==Parkes70)
        return SkyCoverageParkes70(fLDeg, fBDeg);
    if (SrvID==SwinIntermLat)
        return SkyCoverageSwinIL(fLDeg, fBDeg);
    if (SrvID==SwinExtended)
        return SkyCoverageSwinExt(fLDeg, fBDeg);
    if (SrvID==Burgay)
        return SkyCoverageBurgay(fLDeg, fBDeg);
    if (SrvID==GMRT)
        return SkyCoverageGMRT(fLDeg, fBDeg);
    if (SrvID==ParkesMultiBeamALLSKY)
        return SkyCoverageAllSky(fLDeg, fBDeg);
    if (SrvID==ParkesMultiBeamPart)
        return SkyCoveragePartSky(fLDeg, fBDeg);
    if (SrvID==SKA1Mid)
        return SkyCoverageSKA1Mid(fLDeg, fBDeg);
    if (SrvID==SKA1MidLC)
        return SkyCoverageSKA1Mid(fLDeg, fBDeg);
    WARNING("No such survey");
    return false;
}






CSurvey::~CSurvey()
{
	IsCovered=0;
}

#include<iostream>
using namespace std;
// returns DM0
RealType CSurvey::m_DM0()
{
   RealType fWavelength = gCLight/(fFrequency*1e6);
   //cout<< "DDM: "<<1000.0 * fSamplingTime * pow(3e2/fWavelength,3) / (8.3e6 * (fReceiverBandwidth/nChannelsNum ) )<<endl;
   return 1000.0 * fSamplingTime * pow(3e2/fWavelength,3) / (8.3e6 * (fReceiverBandwidth/nChannelsNum ) );
}



//FIXME przepisac potem na bool fun(const Galactic)
/* Geometry checking functions  - for each survey */
bool CSurvey::SkyCoverageParkes70(const RealType fLDeg, const RealType fBDeg)
{

   RealType alpha, delta;
   /* This survey is defined in equatorial coordinates. */
   //gal2eq(fLDeg*M_PI/180, fBDeg*M_PI/180, &alpha, &delta);
   Galactic Gal;
   Gal.fL=fLDeg;
   Gal.fB=fBDeg;
   Equatorial Equ = Gal2Equ(Gal);
   alpha=Equ.fRA;
   delta=Equ.fDC;
   alpha=alpha*12/M_PI;
   delta=delta*180/M_PI;
   // printf("alpha=%.3lf hours, delta=%.3lf deg\n", alpha,delta);

   	if ( delta <=  0   &&   0 <= alpha && alpha <= 24 )
   	{
		return true;
	}
      
	return false;
}


bool CSurvey::SkyCoverageSKA1Mid(const RealType fLDeg, const RealType fBDeg)
{
	return SkyCoverageParkesMB(fLDeg,fBDeg);

}

bool CSurvey::SkyCoverageParkesMB(const RealType fLDeg, const RealType fBDeg) {

    /* if ( fabs(b) <= 5   &&   -100 <= l  && l <= 50 ) */
    if (fabs(fBDeg) <= 5 &&
        ((0 <= fLDeg && fLDeg <= 50) || (260 <= fLDeg && fLDeg < 360))) {
        return true;
    }
    return false;
}


bool CSurvey::SkyCoverageSwinIL(const RealType fLDeg, const RealType fBDeg)
{
   /* if ( 5 <= fabs(b) && fabs(b) <= 15   &&   -100 <= l && l <= 50 ) */
   	if ( 5 <= fabs(fBDeg) && fabs(fBDeg) <= 15 && (  (0 <= fLDeg && fLDeg <= 50) || (260 <= fLDeg && fLDeg < 360) ) )
	{
		return true;
	}
	return false;
}


bool CSurvey::SkyCoverageSwinExt(const RealType fLDeg, const RealType fBDeg)
{
   /* if ( 15 <= fabs(b) && fabs(b) <= 30   &&   -100 <= l && l <= 50 ) */
   if (  ( 15 <= fabs(fBDeg) && fabs(fBDeg) <= 30 && (  (0 <= fLDeg && fLDeg <= 50) || (260 <= fLDeg && fLDeg < 360) ) ))
   {
      return true;
   }
   return false;
}


bool CSurvey::SkyCoverageBurgay(const RealType fLDeg, const RealType fBDeg)
{
   /* if( fabs(b) <= 60   &&   -140 <= l && l <= -100) */
   if( fabs(fBDeg) <= 60  &&  220 <= fLDeg && fLDeg <= 260  )
   {
      return true;
   }
   return false;
}

bool CSurvey::SkyCoverageGMRT(const RealType fLDeg, const RealType fBDeg)
{

   RealType alpha, delta;

   /* This survey is defined in equatorial coordinates. */
   //gal2eq(fLDeg*M_PI/180, fBDeg*M_PI/180, &alpha, &delta);
   Galactic Gal;
   Gal.fL=fLDeg;
   Gal.fB=fBDeg;
   Equatorial Equ = Gal2Equ(Gal);
   alpha=Equ.fRA;
   delta=Equ.fDC;
   alpha=alpha*12/M_PI;
   delta=delta*180/M_PI;
   // printf("alpha=%.3lf hours, delta=%.3lf deg\n", alpha,delta);
 
   if ( delta > -40   &&   0 <= alpha && alpha <= 24 )
   {
	   return true;
   }
   return false;
}

bool CSurvey::SkyCoverageAllSky(const RealType fLDeg, const RealType fBDeg)
{
  	return true;
}

bool CSurvey::SkyCoveragePartSky(const RealType fLDeg, const RealType fBDeg)
{
  	if (fabs(fBDeg) <= 30 && ((0 <= fLDeg && fLDeg <= 90) || (200 <= fLDeg && fLDeg <= 360)) )
	{
		return true;
	}
	return false;
}




