#include<Galaxy.hpp>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//TESTS


void PositionTest()
{
	CGalaxy Gal;
	CAreaDensity AD(-25.,-25.,25.,25.,1024,1024)  ;
	ofstream out("InitPos.dat");
	Coordinates Pos;
	for(int i=0;i<10000;i++)
	{
		*Pos.Cart=Gal.RandomPulsarInitPos();
		out<<Pos.Cart->fX<<" "<<Pos.Cart->fY<<" "<<Pos.Cart->fZ<<endl;
		AD.AddPoint(Pos.Cart->fX, Pos.Cart->fY);
	}
	out.close();
	AD.Densify();
	AD.Write("InitPos.mtx");
	AD.Norm();
	AD.Write("InitPosNorm.mtx");

}

void KickVelocityTest()
{
	CGalaxy Gal;
	ofstream out("KickVelocity.dat");
	Coordinates Pos;
	for(int i=0;i<10000;i++)
	{
		*Pos.Cart=Gal.RandomPulsarInitVelocity();
		out<<Pos.Cart->fX<<" "<<Pos.Cart->fY<<" "<<Pos.Cart->fZ<<endl;
	}
	out.close();
}

#include<iomanip>



void PotentialTest()
{
	ofstream out("GalacticPotentialTest.dat");
	CPulsarDynamics P1;
	P1.m_Position.fX=0.;
	P1.m_Position.fY=0.;
	P1.m_Position.fZ=0.;

	out<<scientific<<setprecision(19);
	for(RealType a(0.);a<=40.;a+=0.01)
	{
		P1.m_Position.fX=a;
		out<<a<<" "<<-P1.GalacticPotential(P1.m_Position)<<endl;
	}
	out.close();
}

void EnergyTest()
{
	CPulsarDynamics P1;
	P1.Init();

	ofstream out("EnergyTest.dat");


	RealType fMaxTime(250.),fDTime(0.0001);
	
	RealType fEKin(0.), fEPot(0.), fETot(0.);
	RealType fMass(1.), fVTot(0.);

	out<<scientific<<setprecision(19);
	for(int i(0);i<=25000;i++)
	{	
		fVTot=Length(P1.m_Velocity);
		fEKin=fMass*fVTot*fVTot*0.5;
		fEPot=fMass*P1.GalacticPotential(P1.m_Position);
		fETot=fEKin - fEPot;
		(P1.m_Position);
		P1.EvolveRK3Adaptive(0.1);
		out<<i*0.1<<" "<<P1.m_Position.fX<<" "<<P1.m_Position.fY<<" "<<P1.m_Position.fZ<<" "<<fEKin<<" "<<-fEPot<<" "<<fETot<<endl;
	}
	out.close();
}

void EvoTest()
{
	CPulsarDynamics P1,P2;
	P1.Init();
	P2.m_Position=P1.m_Position;
	P2.m_Velocity=P1.m_Velocity;

	//ofstream e1("Evo1.dat");
	//ofstream e2("Evo2.dat");

	//e1<<setprecision(19);
	//e2<<setprecision(19);

	RealType fMaxTime(250.),fDTime(0.0001);
	
#if 0
	for(RealType fTime(0.);fTime<=fMaxTime;fTime+=fDTime)
	{
		P1.Evolve(fDTime);
		e1<<P1.m_Position.fX<<" "<<P1.m_Position.fY<<endl;
	}
#endif
	for(int n(0);n<1000000;n++)
	{
		P2.Init();
		for(int i(0);i<=25;i++)
		{	
			P2.EvolveRK3Adaptive(10.);
	//	e2<<P2.m_Position.fX<<" "<<P2.m_Position.fY<<endl;
		}
	//e1.close();
	//e2.close();
	}
}


void MultipleMotionTest()
{
	CGalaxy Gal;
	Cartesian *Pos, *Vel, Acc;
	RealType Rho, Vr, Ar, fLength;
	const int nStars(1000);
	Pos=new Cartesian[nStars];
	Vel=new Cartesian[nStars];
	
	///////////////////////////
	//Init
	for(int i(0);i<nStars;i++)
	{
		Pos[i]=Gal.RandomPulsarInitPos();
		//Cartesian Acc;
		Acc=Gal.GalacticAcceleration(Pos[i]);
		Ar=sqrt(Acc.fX*Acc.fX + Acc.fY*Acc.fY);
		fLength=sqrt(Pos[i].fX*Pos[i].fX + Pos[i].fY*Pos[i].fY);
		//Predkosc Keplerowska
		Vr=sqrt(Ar*fLength);
		Vel[i].fZ=Vel[i].fX=Vel[i].fY=0.;
		Vel[i].fX=-Vr*(-Pos[i].fY/fLength);//orienta troche lepsza ;]
		Vel[i].fY=-Vr*(Pos[i].fX/fLength);
	}

	const RealType fMaxTime(250.), fDTime(0.001);//[Myr]
	int outnum(0),step(0);

	for(RealType fTime(0.);fTime<=fMaxTime;fTime+=fDTime)
	{
		for(int i(0);i<nStars;i++)
		{
			Acc=Gal.GalacticAcceleration(Pos[i]);
			Vel[i] = Vel[i] + Acc*fDTime;
			Pos[i] = Pos[i] + Vel[i]*fDTime;
		}
		//check co 1Myr
		if(step%1000 == 0)
		{
		       	CAreaDensity ADi(-25.,-25.,25.,25.,1024,1024);
			for(int i(0);i<nStars;i++) ADi.AddPoint(Pos[i].fX, Pos[i].fY);
			ADi.Densify();	
			string nazwa="./PosTest/" + std::to_string(static_cast<long long int>(outnum)) + ".mtx";
			ADi.Write(nazwa.c_str());
			outnum++;
		}
		step++;
	}

	delete [] Pos;
	delete [] Vel;
	
}
void SingleMotionTest()
{
	CGalaxy Gal;
	ofstream out("SingleMotion.dat");
	Cartesian Pos;
	Cartesian Vel;
	RealType Rho, Vr, Ar;
	Cartesian Acc;
	Pos.fX=1.;
	Pos.fY=Pos.fZ=0.;
	int nSteps(18);
	RealType fAngle(0.),fDAngle(2.*M_PI/nSteps),fLength;
	out<<setprecision(9);

	///////////////////////////
	//Init
	
	//Cartesian Acc;
	Acc=Gal.GalacticAcceleration(Pos);
	Ar=sqrt(Acc.fX*Acc.fX + Acc.fY*Acc.fY + Acc.fZ*Acc.fZ);
	fLength=sqrt(Pos.fX*Pos.fX + Pos.fY*Pos.fY + Pos.fZ*Pos.fZ);
	//Predkosc Keplerowska
	Vr=sqrt(Ar*fLength);
	Vel.fZ=Vel.fX=Vel.fY=0.;
	Vel.fX=-Vr*(-Pos.fY/fLength);//orienta troche lepsza ;]
	Vel.fY=-Vr*(Pos.fX/fLength);


	//out<<Pos.fX<<" "<<Pos.fY<<" "<<Pos.fZ<<" "<<Ar<<" "<<Acc.fX<<" "<<Acc.fY<<" "<<Acc.fZ<<" "<<Vr<<" "<<Vel.fX<<" "<<Vel.fY<<" "<<Vel.fZ<<endl;
	const RealType fMaxTime(250.), fDTime(0.001), fCheck(0.1);//[Myr]
	int outnum(0),step(0);
	for(RealType fTime(0.);fTime<=fMaxTime;fTime+=fDTime)
	{
		Acc=Gal.GalacticAcceleration(Pos);
		Ar=sqrt(Acc.fX*Acc.fX + Acc.fY*Acc.fY + Acc.fZ*Acc.fZ);
		fLength=sqrt(Pos.fX*Pos.fX + Pos.fY*Pos.fY + Pos.fZ*Pos.fZ);

		Vel = Vel + Acc*fDTime;
		Pos = Pos + Vel*fDTime;
		//out<<Pos.fX<<" "<<Pos.fY<<" "<<Pos.fZ<<" "<<Ar<<" "<<Acc.fX<<" "<<Acc.fY<<" "<<Acc.fZ<<" "<<Vr<<" "<<Vel.fX<<" "<<Vel.fY<<" "<<Vel.fZ<<endl;
			
	}
	out.close();

}


void RotationalCurveTest()
{
	CGalaxy Gal;
	Coordinates Pos;
	Pos.Zero();
	RealType dX=0.01;
	RealType Ar;
	Cartesian Acc;

	ofstream out("RotationalCurve.dat");

	for(Pos.Cart->fX=0.1;Pos.Cart->fX<120.;Pos.Cart->fX+=dX)
	{
		Acc=Gal.GalacticAcceleration(*(Pos.Cart));
		Ar=sqrt(Acc.fX*Acc.fX + Acc.fY*Acc.fY + Acc.fZ*Acc.fZ);
		out<<setprecision(9);
		out<<Pos.Cart->fX<<" "<<1e3*sqrt(Ar*Pos.Cart->fX)<<endl;
	}
	out.close();
	
}

void YusifovKucukRhoDistTest()
{	
	CGalaxy CG;
	RealType fRho(0.);
	const char *nazwa[4]={"Norma.Arm","Carina-Sagittarius.Arm","Perseus.Arm","Crux-Scutum.Arm"};
	for(int i(1);i<=4;i++)
	{	
		fstream out(nazwa[i-1], std::ofstream::out );
		CGalaxy CG;
		RealType fRho(4.), fTheta, fX, fY;
		RealType fDRho(0.001);
		while(fRho<18.)
		{
			fTheta=CG.m_SpiralArmTheta(fRho, i);
			fX=fRho*cos(fTheta);
			fY=fRho*sin(fTheta);
			out<<fRho<<" "<<fTheta<<" "<<fX<<" "<<fY<<endl;
			fRho+=fDRho;
		}
		out.close();
	}
}

int main()
{
	//SpiralArmTest();
	//YusifovKucukRhoDistTest();
	//VelocityTest();
	//PositionTest();
	//SingleMotionTest();
	//MultipleMotionTest();
	//EvoTest();
	//EnergyTest();
	PotentialTest();
	return 0;
}

