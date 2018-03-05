#include<Galaxy.hpp>
#include<PulsarDynamics.hpp>
#include<iomanip>
#include<fstream>
#include<string>
#include<ConfigContainer.hpp>
#include<Statistics.hpp>
#include<MCMC.hpp>

#include<Population.hpp>
using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//TESTS



void InitialPositionTest()
{
	CTestGalaxy Gal;
	Gal.PositionTest();
}



void EnergyTestMultiple()
{

	//ofstream out("EnergyTest.dat");


	RealType fMaxTime(250.),fDTime(0.0001);
	RealType fEKin(0.), fEPot(0.), fETot(0.);
	RealType fMass(1.), fVTot(0.);
	Cartesian L;
	RealType fETotBeg, fLZBeg;
	//out<<scientific<<setprecision(19);
	for(int j(0);j<10000;j++)
	{
		CPulsarDynamics P1;
		P1.Init();
		for(int i(0);i<=25000;i++)
		{
			fVTot=Length(P1.Velocity);
   		    fEKin=fMass*fVTot*fVTot*0.5;
			fEPot=fMass*P1.GalacticPotential(P1.Position);
			fETot=fEKin - fEPot;
			L=1.*CrossProduct(P1.Position,P1.Velocity);//masa 1.

			(P1.Position);
			//P1.m_EvolveRKi3Adaptive(0.1);
			//out<<i*0.001<<" "<<P1.Position.fX<<" "<<P1.Position.fY<<" "<<P1.Position.fZ<<" "<<fEKin<<" "<<-fEPot<<" "<<fETot<<" "<<L.fZ<<endl;
			if(i==0)
			{
				fETotBeg=fETot;
				fLZBeg=L.fZ;
			}


			P1.Evolve(0.1);
		}
		cout<<scientific<<setprecision(19)<<(fETot-fETotBeg)/fETot<<" "<<(L.fZ-fLZBeg)/L.fZ<<endl;
	}
	//out.close();
}


void EvoTest()
{
	CPulsarDynamics P1,P2;
	P1.Init();
	P2.Position=P1.Position;
	P2.Velocity=P1.Velocity;

	//ofstream e1("Evo1.dat");
	//ofstream e2("Evo2.dat");

	//e1<<setprecision(19);
	//e2<<setprecision(19);

	RealType fMaxTime(250.),fDTime(0.0001);

#if 0
	for(RealType fTime(0.);fTime<=fMaxTime;fTime+=fDTime)
	{
		P1.Evolve(fDTime);
		e1<<P1.Position.fX<<" "<<P1.Position.fY<<endl;
	}
#endif
	for(int n(0);n<1000000;n++)
	{
		P2.Init();
		for(int i(0);i<=25;i++)
		{
			//P2.m_EvolveRK3Adaptive(10.);
			P2.Evolve(10.);
	//	e2<<P2.Position.fX<<" "<<P2.Position.fY<<endl;
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

#if 0
int main()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit(); 

	CTestGalaxy *GalaxyTest = new CTestGalaxy;
//	GalaxyTest->RunAll();
	delete GalaxyTest;


	CTestPopulation Pop;
	Pop.WriteReadTest();
//	Pop.Init();
//	Pop.Write2File("Population.Init.mtx");
//	Pop.Evolve();
//	Pop.Write2File("Population.End.mtx");



//	CTestObservations Obs;
//	Obs.Test();


//	CTestModel Model;
//	Model.Test();

//	ArmsTest();
	//SpiralArmTest();
	//YusifovKucukRhoDistTest();
	//VelocityTest();
//	PositionTest();
	//SingleMotionTest();
	//MultipleMotionTest();
	//EvoTest();
	//EnergyTestMultiple();
	//PotentialTest();
	return 0;
}
#endif

int main(int argc, char** argv)
{
	CTestGalaxy Gal;
	Gal.ArmsTest();
	Gal.PositionTest();
	return 0;
	if(argc!=2)
	{
		cout<<"./GalacticTests in_pop.bin"<<endl;
		return 1;
	}
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	pCfgPtr->PostConfigInit();

	CPopulation Pop;
	Pop.ReadBinary(argv[1]);
 	Pop.Write("");

	return 0;
	pCfgPtr->Sim.bDynamics=true;
	pCfgPtr->Sim.bPhysics=true;
	pCfgPtr->Sim.bObservations=true;
	pCfgPtr->Sim.bVerbose=true;    
	pCfgPtr->Sim.nMaxPulsarNum=100000;

	cout<<"in: "<<argv[1]<<endl;
	Pop.Init();
        Pop.Evolve();



	Pop.WriteBinary(argv[1]);
	Pop.Position("Population.XY.mtx");
	Pop.Project("Population.LB.mtx");
	Pop.Write("");

	return 0;
}
