#include<NSNS.hpp>
#include<fstream>
#include<Projections.hpp>
#include<algorithm>


using namespace std;


CNSNSBinaryPopulation::CNSNSBinaryPopulation()
{
	m_pCfgPtr = CConfiguration::GetInstance();
	m_pRanPtr = m_pCfgPtr->pRan;

}

CNSNSBinaryPopulation::~CNSNSBinaryPopulation()
{
	if(m_pBinary!=NULL) delete [] m_pBinary;
}


void CNSNSBinaryPopulation::Read(const char *Filename)
{
	const RealType cfkms2kpcMyr(1e-3);//FIXME poprawic dokladnosc
	ifstream in(Filename);

	//count number of lines
	m_nNumOfBinaries = count(istreambuf_iterator<char>(in), istreambuf_iterator<char>(), '\n');

	in.clear();
	in.seekg(0, ios::beg);
	//create the containter
	m_pBinary = new CBinary[m_nNumOfBinaries];
	//read input
	for(int i(0);i<m_nNumOfBinaries;++i)
	{
		in>>m_pBinary[i].nId;
		in>>m_pBinary[i].fTimeBoring;
		in>>m_pBinary[i].fTime1;
		in>>m_pBinary[i].fTime2;
		in>>m_pBinary[i].Kick1.fX;
		in>>m_pBinary[i].Kick1.fY;
		in>>m_pBinary[i].Kick1.fZ;
		m_pBinary[i].Kick1 = cfkms2kpcMyr*m_pBinary[i].Kick1;
		in>>m_pBinary[i].Kick2.fX;
		in>>m_pBinary[i].Kick2.fY;
		in>>m_pBinary[i].Kick2.fZ;
		m_pBinary[i].Kick2 = cfkms2kpcMyr*m_pBinary[i].Kick2;
	}

	in.close();
}

void CNSNSBinaryPopulation::Write(const char *Filename)
{
	ofstream out(Filename);
	out<<"#Id ";
	out<<"GalX[kpc] ";
	out<<"GalY[kpc] ";
	out<<"GalZ[kpc] ";
	out<<"GalL[deg] ";
	out<<"GalB[deg] ";
	out<<"GalDist[kpc] ";
	out<<"GalVx[kpc/Myr] ";
	out<<"GalVy[kpc/Myr] ";
	out<<"GalVz[kpc/Myr] ";
	out<<endl;



	for(int i(0);i<m_nNumOfBinaries;++i)
	{
		out<<m_pBinary[i].nId<<" ";
		out<<m_pBinary[i].Pos.fX<<" ";
		out<<m_pBinary[i].Pos.fY<<" ";
		out<<m_pBinary[i].Pos.fZ<<" ";
		out<<m_pBinary[i].GalacticPosDeg.fL<<" ";
		out<<m_pBinary[i].GalacticPosDeg.fB<<" ";
		out<<m_pBinary[i].GalacticPosDeg.fDist<<" ";
		out<<m_pBinary[i].Vel.fX<<" ";
		out<<m_pBinary[i].Vel.fY<<" ";
		out<<m_pBinary[i].Vel.fZ<<endl;

	}
	out.close();
}

void CNSNSBinaryPopulation::Init()
{
	//TODO Zmienic w jednolita metoda dla Pulsar + NSNSBinary
	RealType fPositionStartTime(0.);
	for(int i(0);i<m_nNumOfBinaries;++i)
	{
		m_pBinary[i].Pos=RandomPulsarInitPos();
	
		fPositionStartTime = m_pBinary[i].fTimeBoring + m_pBinary[i].fTime1 + m_pBinary[i].fTime2;

		//TODO Sztywna rotacja Galaktyki wstecz - sprawa dyskusyjna!!!
		m_pBinary[i].Pos=RotateVectorXYPlane(m_pBinary[i].Pos, m_ComputeRigidRotationAngle(-fPositionStartTime));

		//Poczatkowa predkosc Keplerowska
		Cartesian Acc;
		Acc=GalacticAcceleration(m_pBinary[i].Pos);
		RealType Vr=sqrt(Length(Acc)*Length(m_pBinary[i].Pos));


		m_pBinary[i].Vel.fZ = m_pBinary[i].Vel.fX = m_pBinary[i].Vel.fY=0.;
		m_pBinary[i].Vel.fX = -Vr*(-m_pBinary[i].Pos.fY/Length(m_pBinary[i].Pos));//orienta troche lepsza ;]
		m_pBinary[i].Vel.fY = -Vr*(m_pBinary[i].Pos.fX/Length(m_pBinary[i].Pos));

	}

}

void CNSNSBinaryPopulation::Evolve()
{
	int nSteps(0);
	RealType fLocalTimeStep;

	for(int i(0);i<m_nNumOfBinaries;++i)
	{

		m_EvolveDynamics(m_pBinary + i, m_pBinary[i].fTimeBoring);
		m_Kick(m_pBinary + i, m_pBinary[i].Kick1);
		m_EvolveDynamics(m_pBinary + i, m_pBinary[i].fTime1);
		m_Kick(m_pBinary + i, m_pBinary[i].Kick2);
		m_EvolveDynamics(m_pBinary + i, m_pBinary[i].fTime2);

	}
}


void CNSNSBinaryPopulation::m_Kick(CBinary *pBinary, const Cartesian Kick) 
{
	pBinary->Vel = pBinary->Vel + Kick;
}

void CNSNSBinaryPopulation::m_EvolveDynamics(CBinary *pBinary, const RealType fEvoTime) 
{
	int nSteps( static_cast<int>(fEvoTime/m_pCfgPtr->Sim.fMaxTimeStep)+1 );
	RealType fLocalTimeStep( fEvoTime/nSteps );

	for(int j(0);j<nSteps;++j) m_EvolvePositionVerlet(pBinary, fLocalTimeStep);
}




void CNSNSBinaryPopulation::m_EvolvePositionVerlet(CBinary *pBinary, const RealType fTimeStep) 
{
	Cartesian Acc;
	
	//first halfstep for the position
	pBinary->Pos  = pBinary->Pos + 0.5*pBinary->Vel*fTimeStep;
	Acc = GalacticAcceleration(pBinary->Pos);
	//velocity full-leap-step
	pBinary->Vel = pBinary->Vel + Acc*fTimeStep;
	//second halfstep for the position
	pBinary->Pos = pBinary->Pos + 0.5*pBinary->Vel*fTimeStep;

}



void CNSNSBinaryPopulation::Project()
{
	for(int i(0);i<m_nNumOfBinaries;++i)
	{
		m_pBinary[i].GalacticPos=Cart2Gal(m_pBinary[i].Pos, m_SunCurrentPosition);
		m_pBinary[i].GalacticPosDeg=m_pBinary[i].GalacticPos;
		m_pBinary[i].GalacticPosDeg.fL *= (180./M_PI);
		m_pBinary[i].GalacticPosDeg.fB *= (180./M_PI);
	}
}

void CNSNSBinaryPopulation::WriteProjectionMtx()
{
	//TODO Check memory constraint
	CProjection Proj(512,256);
        CAreaDensity AD(-25.,-25.,25.,25.,256,256);
	Proj.Zero();
	for(int i(0);i<m_nNumOfBinaries;++i)
	{
		Proj.AddPoint( (m_pBinary[i].GalacticPos.fL>M_PI ?m_pBinary[i].GalacticPos.fL-2*M_PI: m_pBinary[i].GalacticPos.fL) ,m_pBinary[i].GalacticPos.fB);
		AD.AddPoint(m_pBinary[i].Pos.fX, m_pBinary[i].Pos.fY);

	}
	Proj.SphereDensify();
	Proj.Norm();
	Proj.Write("NSNS.Sky.mtx");

        AD.Densify();
        AD.Norm();
        AD.Write("NSNS.GalacticPlane.mtx");

}

