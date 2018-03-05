#include<BasicStar.hpp>

void CBasicStar::PrintDescriptiveStarInfo(std::ostream& out) const
{
	out<<"bBeaming: "<<bBeaming<<std::endl;
	out<<"bInsideGalaxy: "<<bInsideGalaxy<<std::endl;
	out<<"bIsLost: "<<bIsLost<<std::endl;
	out<<"fProbability: "<<fProbability<<std::endl;
	out<<"fAge: "<<fAge<<std::endl;
	out<<"fRelativeAge: "<<fRelativeAge<<std::endl;
	out<<"fMass: "<<fMass<<std::endl;
	out<<"fPeriod: "<<fPeriod<<std::endl;
	out<<"fPeriodDot: "<<fPeriodDot<<std::endl;
	out<<"fBField: "<<fBField<<std::endl;
	out<<"fInitialPeriod: "<<fInitialPeriod<<std::endl;
	out<<"fInitialPeriodDot: "<<fInitialPeriodDot<<std::endl;
	out<<"fInitialBField: "<<fInitialBField<<std::endl;
	out<<"fFinalBField: "<<fFinalBField<<std::endl;
	out<<"fSpectralIndex: "<<fSpectralIndex<<std::endl;
	out<<"fL400: "<<fL400<<std::endl;
	out<<"fS400: "<<fS400<<std::endl;
	out<<"fS1400: "<<fS1400<<std::endl;
	out<<"fDM: "<<fDM<<std::endl;
	out<<"Position.fX: "<<Position.fX<<std::endl;
	out<<"Position.fY: "<<Position.fY<<std::endl;
	out<<"Position.fZ: "<<Position.fZ<<std::endl;
	out<<"Velocity.fX: "<<Velocity.fX<<std::endl;
	out<<"Velocity.fY: "<<Velocity.fY<<std::endl;
	out<<"Velocity.fZ: "<<Velocity.fZ<<std::endl;
	out<<"InitialPosition.fX: "<<InitialPosition.fX<<std::endl;
	out<<"InitialPosition.fY: "<<InitialPosition.fY<<std::endl;
	out<<"InitialPosition.fZ: "<<InitialPosition.fZ<<std::endl;
	out<<"InitialKick.fX: "<<InitialKick.fX<<std::endl;
	out<<"InitialKick.fY: "<<InitialKick.fY<<std::endl;
	out<<"InitialKick.fZ: "<<InitialKick.fZ<<std::endl;
	out<<"GalacticPos.fL: "<<GalacticPos.fL<<std::endl;
	out<<"GalacticPos.fB: "<<GalacticPos.fB<<std::endl;
	out<<"GalacticPos.fDist: "<<GalacticPos.fDist<<std::endl;
	out<<"GalacticPosDeg.fL: "<<GalacticPosDeg.fL<<std::endl;
	out<<"GalacticPosDeg.fB: "<<GalacticPosDeg.fB<<std::endl;
	out<<"GalacticPosDeg.fDist: "<<GalacticPosDeg.fDist<<std::endl;
}


void CBasicStar::PrintStarInfo(std::ostream& out) const
{
	out<<bBeaming<<" ";
	out<<bInsideGalaxy<<" ";
	out<<bIsLost<<" ";
	out<<fProbability<<" ";
	out<<fAge<<" ";
	out<<fRelativeAge<<" ";
	out<<fMass<<" ";
	out<<fPeriod<<" ";
	out<<fPeriodDot<<" ";
	out<<fBField<<" ";
	out<<fInitialPeriod<<" ";
	out<<fInitialPeriodDot<<" ";
	out<<fInitialBField<<" ";
	out<<fFinalBField<<" ";
	out<<fSpectralIndex<<" ";
	out<<fL400<<" ";
	out<<fS400<<" ";
	out<<fS1400<<" ";
	out<<fDM<<" ";
	out<<Position.fX<<" ";
	out<<Position.fY<<" ";
	out<<Position.fZ<<" ";
	out<<Velocity.fX<<" ";
	out<<Velocity.fY<<" ";
	out<<Velocity.fZ<<" ";
	out<<InitialPosition.fX<<" ";
	out<<InitialPosition.fY<<" ";
	out<<InitialPosition.fZ<<" ";
	out<<InitialKick.fX<<" ";
	out<<InitialKick.fY<<" ";
	out<<InitialKick.fZ<<" ";
	out<<GalacticPos.fL<<" ";
	out<<GalacticPos.fB<<" ";
	out<<GalacticPos.fDist<<" ";
	out<<GalacticPosDeg.fL<<" ";
	out<<GalacticPosDeg.fB<<" ";
	out<<GalacticPosDeg.fDist<<" ";
}


void CBasicStar::ResetToInit(const bool bPhysics, const bool bDynamics)
{

	fRelativeAge=-fAge;
	fAge=fMass=0.;
	if(bPhysics)
	{
		bBeaming=false;//random
		fProbability=0.;
		fPeriod=fInitialPeriod;
		fPeriodDot=fInitialPeriodDot;
		fBField=fInitialBField;
		fSpectralIndex=0.;
		fL400=fS400=fS1400=0.;
	}
	if(bDynamics)
	{
		bInsideGalaxy=true;
		bIsLost=false;
		Position=InitialPosition;
		//Velocity= //with kick it's random
	}
}
