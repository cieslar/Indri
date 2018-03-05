#include<PulsarDynamics.hpp>
#include<ConfigContainer.hpp>
#include<limits>



Cartesian CPulsarDynamics::GetPosition() const
{
	Cartesian Pos=Position;
	return Pos;
}


void CPulsarDynamics::Init(const RealType fPositionStartTime)
{
	CConfiguration* pCfgPtr = CConfiguration::GetInstance();
	fAge=0.;
//init
	//WYlosuj poczatkowa pozycje	
	Position=RandomPulsarInitPos();
	Position=RotateVectorXYPlane(Position, m_ComputeRigidRotationAngle(fPositionStartTime));
	
	InitialPosition=Position;

	//Poczatkowa predkosc Keplerowska
	Cartesian Acc;
	Acc=GalacticAcceleration(Position);
	RealType Vr=sqrt(Length(Acc)*Length(Position));

	Velocity.fZ = Velocity.fX=Velocity.fY=0.;
	Velocity.fX = -Vr*(-Position.fY/Length(Position));//orienta troche lepsza ;]
	Velocity.fY = -Vr*(Position.fX/Length(Position));
	//No i jeszcze kick
	if(pCfgPtr->Sim.KickingMethod == PaczynskiKicking)
	{
		InitialKick = PaczynskiKick();
	}
	if(pCfgPtr->Sim.KickingMethod == HobbsKicking)
	{
		InitialKick = HobbsKick();	
	}
	if(pCfgPtr->Sim.KickingMethod != HobbsKicking && pCfgPtr->Sim.KickingMethod != PaczynskiKicking) ERROR("Nie chce mi sie pisac funkcji-switcha, ze ja pierdole...");
	Velocity = Velocity + InitialKick;

}

#if 0
void CPulsarDynamics::ProjectTest()
{
	Cartesian Pos;
	Galactic Gal;
	Pos.fZ=0.;
	Pos.fY=0.;
	for(double a=-25;a<=25;a+=0.1)
	{
		Pos.fX=a;
		Gal=Cart2Gal(Pos,m_SunCurrentPosition);
		cout<<Pos.fX<<" "<<Gal.fL<<" "<<Gal.fB<<" "<<Gal.fDist<<endl;

	}
}
#endif

Galactic CPulsarDynamics::m_Project(const RealType fProjectionTime)
{
	//return Cart2Gal(Position, m_SunCurrentPosition, m_ComputeRigidRotationAngle(fProjectionTime) );
	//return Cart2Gal(Position, Pos, 0. );
	return Cart2Gal(Position,m_SunCurrentPosition);
}

#if 0
//to jest wylcznie dla starego leapfroga
void CPulsarDynamics::EvolveInit(const RealType fMaxTimeStep)
{
	m_VelocityHalfLeap(fMaxTimeStep);
}

void CPulsarDynamics::EvolveEnd(const RealType fMaxTimeStep)
{
	m_VelocityHalfLeap(-fMaxTimeStep);
}
#endif

void CPulsarDynamics::Evolve(const RealType fMaxTimeStep)
{
	//m_EvolveRK3Adaptive(fMaxTimeStep);
	//m_EvolveNewtonForward(fMaxTimeStep);
	//m_EvolveLeapfrog(fMaxTimeStep);
	m_EvolvePositionVerlet(fMaxTimeStep);
	//m_EvolvePositionVerletDynamic(fMaxTimeStep);
}

void CPulsarDynamics::m_EvolveNewtonForward(const RealType fTimeStep)
{
	Cartesian Acc=GalacticAcceleration(Position);

	Velocity = Velocity + Acc*fTimeStep;
	Position = Position + Velocity*fTimeStep;
}

void CPulsarDynamics::m_VelocityHalfLeap(const RealType fTimeStep)
{
	Cartesian Acc=GalacticAcceleration(Position);
	Velocity = Velocity + Acc*(fTimeStep/2.);

}

void CPulsarDynamics::m_EvolveLeapfrog(const RealType fTimeStep)
{

	Cartesian Acc;
	
	Position = Position + Velocity*fTimeStep;
	Acc=GalacticAcceleration(Position);
	Velocity = Velocity + Acc*fTimeStep;
}

void CPulsarDynamics::m_EvolvePositionVerlet(const RealType fTimeStep)
{
	Cartesian Acc;
	
	//first halfstep for the position
	Position = Position + 0.5*Velocity*fTimeStep;
	Acc = GalacticAcceleration(Position);
	//velocity full-leap-step
	Velocity = Velocity + Acc*fTimeStep;
	//second halfstep for the position
	Position = Position + 0.5*Velocity*fTimeStep;



}

void CPulsarDynamics::m_EvolvePositionVerletDynamic(const RealType fMaxTimeStep)
{
	Cartesian Acc;
	Cartesian PosBeg, PosEnd, VelBeg, VelEnd;

	PosBeg = Position;
	VelBeg = Velocity;
	RealType fTimeStep = fMaxTimeStep;
	
	//first halfstep for the position
	PosEnd = PosBeg + 0.5*VelBeg*fTimeStep;
	Acc = GalacticAcceleration(PosEnd);
	//velocity full-leap-step
	VelEnd = VelBeg + Acc*fTimeStep;
	//second halfstep for the position
	PosEnd = PosEnd + 0.5*VelEnd*fTimeStep;

	PosEnd = PosEnd - PosBeg;

	if(Length(PosEnd)/(Length(VelEnd)*fTimeStep) < 5. )
	{
		int n = static_cast<int>( 5.*Length(VelEnd)*fMaxTimeStep/Length(PosEnd) ) + 1;
		fTimeStep = fMaxTimeStep/n;
		
		PosBeg = Position;
		VelBeg = Velocity;

		
		
		for(int i(0);i<n;i++)
		{
			PosEnd = PosBeg + 0.5*VelBeg*fTimeStep;
			Acc = GalacticAcceleration(PosEnd);
			//velocity full-leap-step
			VelEnd = VelBeg + Acc*fTimeStep;
			//second halfstep for the position
			PosEnd = PosEnd + 0.5*VelEnd*fTimeStep;
			
			VelBeg=VelEnd;
			PosBeg=PosEnd;
		}
	}
	

	Position=PosEnd;
	Velocity=VelEnd;




}




//Nie zrozumialem do konca ktorego to jest rzedu metoda
//http://www.artcompsci.org/kali/vol/runge_kutta/ch05.html
void CPulsarDynamics::m_EvolveRK3Adaptive(const RealType fMaxTimeStep)
{
	RealType fTimeStep = fMaxTimeStep;
	RealType fAbsoluteEps     =  1e6*std::numeric_limits<RealType>::epsilon();
	RealType fError;
	Cartesian Acc;
	Cartesian PreviousPos, PreviousVel;
	Cartesian Pos, Vel, TempPos, K1, K2, K3;
	Cartesian DPos;

	Pos=Position;
	Vel=Velocity;
		

	K1=GalacticAcceleration(Pos);

	TempPos=Pos + 0.5*Vel*fTimeStep + K1*fTimeStep*fTimeStep*(1./8.);
	K2=GalacticAcceleration(TempPos);
		
	TempPos=Pos + Vel*fTimeStep + 0.5*K2*fTimeStep*fTimeStep;
	K3=GalacticAcceleration(TempPos);

	Pos = Pos + Vel*fTimeStep + (K1 + 2.*K2)*fTimeStep*fTimeStep*(1./6.);
	Vel = Vel + (K1 + 4.*K2 + K3)*fTimeStep*(1./6.);

	int i(1);
	do
	{
		i<<=1;
		fTimeStep/=2.;

		PreviousPos=Pos;
		PreviousVel=Vel;
	
		Pos=Position;
		Vel=Velocity;
		for(int j(0);j<i;j++)
		{
			K1=GalacticAcceleration(Pos);

			TempPos=Pos + 0.5*Vel*fTimeStep + K1*fTimeStep*fTimeStep*(1./8.);
			K2=GalacticAcceleration(TempPos);
		
			TempPos=Pos + Vel*fTimeStep + 0.5*K2*fTimeStep*fTimeStep;
			K3=GalacticAcceleration(TempPos);

			Pos = Pos + Vel*fTimeStep + (K1 + 2.*K2)*fTimeStep*fTimeStep*(1./6.);
			Vel = Vel + (K1 + 4.*K2 + K3)*fTimeStep*(1./6.);

			//FIXME dopisac poprawne sprawdzanie dokladnosci x i v

			

		}
		DPos=Pos - PreviousPos;
		fError = pow(fTimeStep,4.);
	}
	while( !( Length(DPos) < fAbsoluteEps/Length(Pos) ) && !( Length(DPos) < fError/Length(Pos) ) );


	Position=Pos;
	Velocity=Vel;
}




