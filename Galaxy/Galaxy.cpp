#include<Galaxy.hpp>
#include<ConfigContainer.hpp>
#include<cmath>
#include<Debug.hpp>

#include<iomanip>
#include<fstream>
#include<string>

#include<Statistics.hpp>
CGalaxy::CGalaxy() : m_SunCurrentPosition(0., 8.5, 0.)
{
	CConfiguration *pCfgPtr=CConfiguration::GetInstance();
	m_pRanPtr=pCfgPtr->pRan;

	//m_SunCurrentPosition.fX=0.;
	//m_SunCurrentPosition.fY=8.5;//[kpc]
	//m_SunCurrentPosition.fZ=0.;
	m_SunOmega=m_ComputeSunOmega();//zakladam sztywna rotacje

}

CGalaxy::~CGalaxy()
{

}

inline void CGalaxy::m_ComputePartialElements(const Cartesian Position) 
{
	fRho = sqrt (Position.fX*Position.fX+Position.fY*Position.fY);
	fR=sqrt(Position.fX*Position.fX + Position.fY*Position.fY + Position.fZ*Position.fZ);
//	fR2 = fR*fR;


}
///Miyamoto & Nagai potential for Disk 
///Mass[1e10Msun] 
///fParA,fParB,Position[kpc]
inline RealType CGalaxy::m_MiyamotoNagaiDisk(const RealType fMass, const RealType fParA, const RealType fParB, const Cartesian Position) const
{
//	RealType fRho = sqrt (Position.fX*Position.fX+Position.fY*Position.fY);
	return gfG*fMass/sqrt( fRho*fRho + pow( fParA + sqrt(fParB*fParB+Position.fZ*Position.fZ) , 2. ));
}


inline Cartesian CGalaxy::m_DivMiyamotoNagaiDisk(const RealType fMass, const RealType fParA, const RealType fParB, const Cartesian Position) const
{
//	RealType fRho=sqrt(Position.fX*Position.fX + Position.fY*Position.fY);

	RealType CommonPart=gfG*fMass/pow( fRho*fRho + pow( fParA + sqrt(fParB*fParB+Position.fZ*Position.fZ) , 2. ),3./2.);


	Cartesian Div(Position);
//	Div.fX=Div.fY=Div.fZ=0.;
//	Div.fX=-CommonPart*Position.fX;
//	Div.fY=-CommonPart*Position.fY;
//	Div.fZ=-CommonPart*Position.fZ*(1.+fParA/sqrt(Position.fZ*Position.fZ+fParB*fParB));
	Div.fZ *= (1.+fParA/sqrt(Position.fZ*Position.fZ+fParB*fParB));
	return -CommonPart*Div;
}


inline RealType  CGalaxy::m_MiyamotoNagaiBuldge(const RealType fMass, const RealType fParB, const Cartesian Position) const
{
//	RealType fR=sqrt(Position.fX*Position.fX + Position.fY*Position.fY + Position.fZ*Position.fZ);
	return gfG*fMass/sqrt( fR*fR +  fParB*fParB);
}


inline Cartesian CGalaxy::m_DivMiyamotoNagaiBuldge(const RealType fMass, const RealType fParB, const Cartesian Position) const
{
//	RealType fR=sqrt(Position.fX*Position.fX + Position.fY*Position.fY + Position.fZ*Position.fZ);
	RealType CommonPart=gfG*fMass/pow( fR*fR +  fParB*fParB, 3./2.);
	
	Cartesian Div(Position);
//	Div.fX=Div.fY=Div.fZ=0.;
//	Div.fX=-CommonPart*Position.fX;
//	Div.fY=-CommonPart*Position.fY;
//	Div.fZ=-CommonPart*Position.fZ;
	return -CommonPart*Div;
}


//Potential of the Hernquist's Shperoid
//fMass[1e10Msun]
//fParC,Position[kpc]
inline RealType CGalaxy::m_Hernquist(const RealType fMass, const RealType fParC, const Cartesian Position) const
{
//	RealType fR=sqrt(Position.fX*Position.fX + Position.fY*Position.fY + Position.fZ*Position.fZ);
	return gfG*fMass/(fR+fParC);
}

inline Cartesian CGalaxy::m_DivHernquist(const RealType fMass, const RealType fParC, const Cartesian Position) const
{
//	RealType fR=sqrt(Position.fX*Position.fX + Position.fY*Position.fY + Position.fZ*Position.fZ);
	RealType CommonPart=gfG*fMass/(fR*(fR+fParC)*(fR+fParC));
	
	Cartesian Div(Position);
	//Div.fX=Div.fY=Div.fZ=0.;
	//Div.fX=-CommonPart*Position.fX;
	//Div.fY=-CommonPart*Position.fY;
	//Div.fZ=-CommonPart*Position.fZ;
	return -CommonPart*Div;
}

//Potential of Paczynski's Halo
//fMass[1e10Msun]
//
inline RealType CGalaxy::m_Paczynski(const RealType fMass, const RealType fRCore, const RealType fRCut, const Cartesian Position) const
{
//	RealType fR=sqrt(Position.fX*Position.fX + Position.fY*Position.fY + Position.fZ*Position.fZ);
	if(fR<fRCut)
	{
		return -(gfG*fMass/fRCore)*(0.5*log(1.+ fR*fR/(fRCore*fRCore)) + (fRCore/fR)*atan(fR/fRCore));
	}
	else
	{
		return -gfG*fMass/fR + gfG*fMass/fRCut -(gfG*fMass/fRCore)*(0.5*log(1.+ fRCut*fRCut/(fRCore*fRCore)) + (fRCore/fRCut)*atan(fRCut/fRCore));
	}
}

inline Cartesian CGalaxy::m_DivPaczynski(const RealType fMass, const RealType fRCore, const RealType fRCut, const Cartesian Position) const
{

//	RealType fR=sqrt(Position.fX*Position.fX + Position.fY*Position.fY + Position.fZ*Position.fZ);
	RealType CommonPart;
	if(fR<fRCut)
	{
		CommonPart=(gfG*fMass/(fR*fR*fR))*(atan(fR/fRCore) - fR/fRCore);
	}
	else
	{
		CommonPart=gfG*fMass/(fR*fR*fR);
	}
	
	Cartesian Div(Position);
	//Div.fX=CommonPart*Position.fX;
	//Div.fY=CommonPart*Position.fY;
	//Div.fZ=CommonPart*Position.fZ;
	return CommonPart*Div;
}




RealType CGalaxy::GalacticPotential(const Cartesian Position) 
{
//FIXME przeniesc parametry do configu
	//Paczynski parameters
	RealType fMassHalo=5.;
	RealType fRCut=100.;
	RealType fRCore=6.;

	//MyamotoNagai parameters for Buldge
	RealType fMassBuldge=1.12;
	RealType fParBBuldge=0.277;
	//for Disk
	RealType fMassDisk=8.78;
	RealType fParADisk=4.2;
	RealType fParBDisk=0.198;

	m_ComputePartialElements(Position);

	RealType Potential= m_Paczynski(fMassHalo, fRCore, fRCut, Position)
	                  + m_MiyamotoNagaiBuldge(fMassBuldge,fParBBuldge, Position)
	                  + m_MiyamotoNagaiDisk(fMassDisk, fParADisk, fParBDisk, Position);

	return Potential;
}

Cartesian CGalaxy::GalacticAcceleration(const Cartesian Position) 
{
	//Paczynski parameters
	RealType fMassHalo=5.;
	RealType fRCut=100.;
	RealType fRCore=6.;

	//MyamotoNagai parameters for Buldge
	RealType fMassBuldge=1.12;
	RealType fParBBuldge=0.277;
	//for Disk
	RealType fMassDisk=8.78;
	RealType fParADisk=4.2;
	RealType fParBDisk=0.198;

	m_ComputePartialElements(Position);

	Cartesian Acc(   m_DivPaczynski(fMassHalo, fRCore, fRCut, Position)
		       + m_DivMiyamotoNagaiBuldge(fMassBuldge, fParBBuldge, Position)
		       + m_DivMiyamotoNagaiDisk(fMassDisk, fParADisk, fParBDisk, Position) 
		     );
	return Acc;
}



RealType CGalaxy::m_YusifovKucukRhoDist(const RealType fRho) const
{	
	RealType fNorm=1.;//do przecałkowania
	RealType fParamA=1.64;
	RealType fParamB=4.01;
	RealType fParamRho1=0.55;
	RealType fParamRhoSun=8.5;
	return fNorm * pow( (fRho + fParamRho1) / (fParamRhoSun + fParamRho1) , fParamA ) * exp( -fParamB * (fRho-fParamRhoSun)/(fParamRhoSun + fParamRho1) )*2.*M_PI*fRho;
}

RealType CGalaxy::m_RandomPulsarRho() const
{
	RealType fRho,fProbab;
	RealType fProbabMax=m_YusifovKucukRhoDist(5.62863);//z wykresu
	RealType fRhoMax=20.;//[kpc]
	do
	{
		fRho=m_pRanPtr->Flat()*fRhoMax;
		fProbab=m_pRanPtr->Flat()*fProbabMax;
	} 
	while(fProbab>m_YusifovKucukRhoDist(fRho));
	return fRho;
}

RealType CGalaxy::m_SpiralArmTheta(const RealType fRho,
                                   const int nArmNum) const {
    RealType fK(0.), fRho0(0.), fTheta0(0.);
    RealType fTheta(0.);

    // z pracy Faucher-Giguere
    if (nArmNum == NormaArm) {
        fK = 4.25;
        fRho0 = 3.48;
        fTheta0 = 1.57;
        return fK * log(fRho / fRho0) + fTheta0;
    }

    if (nArmNum == CarinaSagittariusArm) {
        fK = 4.25;
        fRho0 = 3.48;
        fTheta0 = 4.71;
        return fK * log(fRho / fRho0) + fTheta0;
    }

    if (nArmNum == PerseusArm) {
        fK = 4.89;
        fRho0 = 4.90;
        fTheta0 = 4.09;
        return fK * log(fRho / fRho0) + fTheta0;
    }

    if (nArmNum == CruxScutumArm) {
        fK = 4.89;
        fRho0 = 4.90;
        fTheta0 = 0.95;
        return fK * log(fRho / fRho0) + fTheta0;
    }

    WARNING("No such arm.");
    return std::numeric_limits<RealType>::signaling_NaN();		
}



//FIXME nazwa -> RandomArmPos
Cartesian CGalaxy::RandomPulsarInitPos() const
{
	Cylindrical Cyli;
	Cartesian Cart;
	Cart.fX=Cart.fY=Cart.fZ=0.;
	//losowanie ramienia
	int nArmNum;
	do 
	{
		nArmNum=static_cast<int>(m_pRanPtr->Flat()*4.+1.);//FIXME - funkcja losujaca dla intow z odpowiednim przedzialem
	}
	while (nArmNum == 5.);

	//losowanie Rho
	Cyli.fRho=m_RandomPulsarRho();

	//obliczenie Theta
	Cyli.fTheta=m_SpiralArmTheta(Cyli.fRho, nArmNum);

	//Zaburzenie Theta
	RealType fThetaCorr = m_pRanPtr->Flat()*2.*M_PI;//FIXME - funkcja losujaca z otwartym przedzialem gornym
	Cyli.fTheta += fThetaCorr * exp(-0.35*Cyli.fRho);

	Cart=Cyli2Cart(Cyli);


#define ROZMYCIE_ORIG
//#define ROZMYCIE_GAUSS_2D

#ifdef ROZMYCIE_ORIG
	RealType fTranslationRho = m_pRanPtr->Gauss(0.07*Cyli.fRho,0.);
	RealType fTranslationTheta = m_pRanPtr->Flat()*M_PI;//chyba potrzebeuje jedynie M_PI skoro losuje tez ujemne Rho
	
	Cart.fX += fTranslationRho*cos(fTranslationTheta);
	Cart.fY += fTranslationRho*sin(fTranslationTheta);
#endif


#ifdef ROZMYCIE_GAUSS_2D
	Cartesian Rozmycie = m_pRanPtr->Gauss2D(0.07*Cyli.fRho,0.);
	Cart = Cart + Rozmycie;
#endif
	//Wylosowanie Z
	RealType fZ0=0.05;//[kpc]
	Cart.fZ=(((m_pRanPtr->Flat()*2. -1.) < 0.) ? -1. : 1.) * m_pRanPtr->Exponential(1./fZ0);//FIXME - dopisac funkcje sgn
	return Cart;
}

RealType CGalaxy::m_PaczysnkiKicksDist(const RealType fVelocity) const
{
	const RealType cfVelParam(560.);
	return 4./( M_PI * cfVelParam * pow(1.+ pow(fVelocity/cfVelParam, 2.), 2.)  );///< [km/s]
}

const RealType gcfkms2kpcMyr(1e-3);//FIXME poprawic dokladnosc

Cartesian CGalaxy::HobbsKick() const
{
	const RealType fAParam(265.);//km/s
	return m_pRanPtr->UniformShpereDirection()*m_pRanPtr->Maxwell(fAParam)*gcfkms2kpcMyr;
}



Cartesian CGalaxy::PaczynskiKick() const
{
	
	RealType fVel,fProbab;
	RealType fProbabMax=m_PaczysnkiKicksDist(0.);
	RealType fVelMax=1400.;//[km/s]

	do
	{
		fVel=m_pRanPtr->Flat()*fVelMax;
		fProbab=m_pRanPtr->Flat()*fProbabMax;
	} 
	while(fProbab>m_PaczysnkiKicksDist(fVel));
	//const RealType cfkms2kpcMyr(1e-3);//FIXME poprawic dokladnosc
	fVel*=gcfkms2kpcMyr;//w [kpc/Myr]
	return m_pRanPtr->UniformShpereDirection()*fVel;
}

RealType CGalaxy::m_ComputeSunOmega() 
{
	Cartesian SunVelocity;
	Cartesian Acc;
	Acc=GalacticAcceleration(m_SunCurrentPosition);

	//Predkosc Keplerowska
	RealType Vr=sqrt(Length(Acc)*Length(m_SunCurrentPosition));

	//m_Velocity.fZ = Vel.fX=Vel.fY=0.;//w konstruktorze
	SunVelocity.fX = -Vr*(-m_SunCurrentPosition.fY/Length(m_SunCurrentPosition));//orienta troche lepsza ;]
	SunVelocity.fY = -Vr*(m_SunCurrentPosition.fX/Length(m_SunCurrentPosition));
	SunVelocity.fZ = 0.;

	return Length(SunVelocity)/Length(m_SunCurrentPosition);
}

RealType CGalaxy::m_ComputeRigidRotationAngle(const RealType fTimePoint) const
{
	return fTimePoint*m_SunOmega;
}











/////////////////////////TESTS//////////////////////////////

CTestGalaxy::CTestGalaxy() 
{

}

void CTestGalaxy::RunAll()
{
	ArmsTest();
	KickVelocityTest();
	PositionTest();
	CoordinatesTransformTest();
	PotentialTestDirection();
	PotentialTest();



}





///TODO tutaj dopisze tez test gestosci ramion
void CTestGalaxy::ArmsTest()
{	
	cout<<"Arms test"<<endl;	
	Cartesian Cart;
	Cylindrical Cyli;

	const RealType fDRho(0.001);
	const char *nazwa[4]={"TestGalaxy.NormaArm.dat","TestGalaxy.CarinaSagittariusArm.dat","TestGalaxy.PerseusArm.dat","TestGalaxy.CruxScutumArm.dat"};
	for(int i(1);i<=4;i++)
	{	
		fstream out(nazwa[i-1], std::ofstream::out );
		out<<"# Rho[kpc] Theta[rad] Z[kpc] X[kpc] Y[kpc] Z[kpc]"<<endl;

		Cyli.fRho=3.;
		while(Cyli.fRho<30.)
		{
			Cyli.fTheta=m_SpiralArmTheta(Cyli.fRho, i);
			Cart=Cyli2Cart(Cyli);
			out<<Cyli.fRho<<" "<<Cyli.fTheta<<" "<<Cyli.fZ<<" "<<Cart.fX<<" "<<Cart.fY<<" "<<Cart.fZ<<endl;
			Cyli.fRho+=fDRho;
		}
		out.close();
	}
}

void CTestGalaxy::PaczynksiKickTest(const int nNum)
{
	ofstream out("TestGalaxy.PaczynskiKick.dat");
	Cartesian Pos;
	out<<"# Vx[kpc/Myr] Vy[kpc/Myr] Vz[kpc/Myr] Vtot[kpc/Myr]"<<endl;
	for(int i=0;i<nNum;i++)
	{
		Pos=PaczynskiKick();
		out<<Pos.fX<<" "<<Pos.fY<<" "<<Pos.fZ<<" "<<sqrt(Pos.fX*Pos.fX + Pos.fY*Pos.fY + Pos.fZ*Pos.fZ)<<endl;
	}
	out.close();
}


void CTestGalaxy::HobbsKickTest(const int nNum)
{
	ofstream out("TestGalaxy.HobbsKick.dat");
	Cartesian Pos;
	out<<"# Vx[kpc/Myr] Vy[kpc/Myr] Vz[kpc/Myr] Vtot[kpc/Myr]"<<endl;
	for(int i=0;i<nNum;i++)
	{
		Pos=HobbsKick();
		out<<Pos.fX<<" "<<Pos.fY<<" "<<Pos.fZ<<" "<<sqrt(Pos.fX*Pos.fX + Pos.fY*Pos.fY + Pos.fZ*Pos.fZ)<<endl;
	}
	out.close();
}





void CTestGalaxy::KickVelocityTest()
{
	cerr<<"Kick velocity test"<<endl;
	ofstream out("TestGalaxy.KickVelocity.dat");
	Cartesian Pos;
	out<<"# Vx[kpc/Myr] Vy[kpc/Myr] Vz[kpc/Myr]"<<endl;
	for(int i=0;i<10000;i++)
	{
		Pos=HobbsKick();
		out<<Pos.fX<<" "<<Pos.fY<<" "<<Pos.fZ<<endl;
	}
	out.close();
}

//FIXME do rozbudowy wraz z nowymi pojemnikami na gestosc
void CTestGalaxy::PositionTest()
{
	cout<<"Initial position test"<<endl;

	CAreaDensity AD(-25.,-25.,25.,25.,512,512);
	//CCubicDensity CD(-25.,-25.,-5.,25.,25.,5.,256,256,256);
	ofstream out("TestGalaxy.InitPos.dat");
	Cartesian Pos;
	const int nNumTested(10000000);

	out<<"# X[kpc] Y[kpc] Z[kpc]"<<endl;
	for(int i=0;i<nNumTested;i++)
	{
		Pos=RandomPulsarInitPos();
		out<<Pos.fX<<" "<<Pos.fY<<" "<<Pos.fZ<<endl;
		AD.AddPoint(Pos.fX, Pos.fY);
		//CD.AddPoint(Pos.fX, Pos.fY,Pos.fZ);
	}
	out.close();
	//AD.Densify();
	AD.Write("TestGalaxy.InitPos.10M.mtx");
//	AD.Norm();
//	AD.Write("TestGalaxy.InitPosNorm.mtx");

//	CD.Densify();
	//CD.Write("TestGalaxy.InitPos3D.10M.dat");

}

//TODO
void CTestGalaxy::CoordinatesTransformTest()
{


}

void CTestGalaxy::PotentialTestDirection()
{
	cout<<"Directional potential test"<<endl;
	ofstream out("TestGalaxy.PotentialTestDirection.dat");

	Cartesian Acc, Pos={0.,0.,0.};
	Cartesian dPos={0.01,0.,0.};//do przebadania kierunku

	RealType Ar,R;
	out<<"# X[kpc] Y[kpc] Z[kpc] -Phi[kpc^2/Myr^2] acc[kpc/Myr^2] Vkep[km/s]"<<endl;
	out<<scientific<<setprecision(19);
	do	
	{
		R   = sqrt(Pos.fX*Pos.fX + Pos.fY*Pos.fY + Pos.fZ*Pos.fZ);
		Acc = GalacticAcceleration(Pos);
		Ar  = sqrt(Acc.fX*Acc.fX + Acc.fY*Acc.fY + Acc.fZ*Acc.fZ);

		out<<Pos.fX<<" "<<Pos.fY<<" "<<Pos.fZ<<" "<<-GalacticPotential(Pos)<<" "<<Length(Acc)<<" "<<1e3*sqrt(Ar*R)<<endl;

		Pos=Pos+dPos;

	}
	while( R <= 30.);
	out.close();
}

void CTestGalaxy::PotentialTest()
{
	cout<<"Potential components test"<<endl;
	ofstream out("TestGalaxy.PotentialTest.dat");
	
	Cartesian Acc, Pos={0.,0.,0.};
	Cartesian dPos={0.01,0.,0.};

	Cartesian DivHalo, DivBuldge, DivDisk, DivTotal;
	RealType Halo, Buldge, Disk, Total;

	//parametry wyciagniete z funkcji 
	//Paczynski parameters
	RealType fMassHalo=5.;
	RealType fRCut=100.;
	RealType fRCore=6.;

	//MyamotoNagai parameters for Buldge
	RealType fMassBuldge=1.12;
	RealType fParBBuldge=0.277;
	//for Disk
	RealType fMassDisk=8.78;
	RealType fParADisk=4.2;
	RealType fParBDisk=0.198;


	out<<"# X[kpc] -PhiDisk[kpc^2/Myr^2] -PhiBuldge -PhiHalo -Phi"<<endl; 
	out<<scientific<<setprecision(19);
	do
	{

		DivHalo   =  m_DivPaczynski(fMassHalo, fRCore, fRCut, Pos);
		DivBuldge =  m_DivMiyamotoNagaiBuldge(fMassBuldge, fParBBuldge, Pos);
		DivDisk   =  m_DivMiyamotoNagaiDisk(fMassDisk, fParADisk, fParBDisk, Pos);
		DivTotal = GalacticAcceleration(Pos);


		Disk   = m_MiyamotoNagaiDisk(fMassDisk, fParADisk, fParBDisk, Pos);
		Buldge = m_MiyamotoNagaiBuldge(fMassBuldge,fParBBuldge, Pos);
		Halo   = m_Paczynski(fMassHalo, fRCore, fRCut, Pos);
		Total  = GalacticPotential(Pos);


		out<<Pos.fX<<" "<<-Disk<<" "<<-Buldge<<" "<<-Halo<<" "<<-Total<<" ";
		out<<Length(DivDisk)<<" "<<Length(DivBuldge)<<" "<<Length(DivHalo)<<" "<<Length(DivTotal)<<" ";
		out<<1e3*sqrt(Length(Pos)*Length(DivDisk))<<" "<<1e3*sqrt(Length(Pos)*Length(DivBuldge))<<" "<<1e3*sqrt(Length(Pos)*Length(DivHalo))<<" "<<1e3*sqrt(Length(Pos)*Length(DivTotal))<<endl;




		Pos=Pos+dPos;

	}
	while( Length(Pos) <= 30.);
	out.close();

}
