#include<Population.hpp>
#include<fstream>
#include<iomanip>
#include<MyMath.hpp>
#include <H5Cpp.h>
#include<sstream>



void CPopulation::WriteTxtDM(const char *FileName) const
{
	ofstream out(FileName);
	for(int i(0);i<m_nNum;i++)
	{
		out<<i<<" "<<pPulsar[i].GalacticPosDeg.fL<<" "<<pPulsar[i].GalacticPosDeg.fB<<" "<<pPulsar[i].GalacticPosDeg.fDist<<" "<<pPulsar[i].fDM<<endl;
	}
	out.close();
}

void CPopulation::ReadTxtDM(const char *FileName)
{
	cerr<<"Warning! Reading YMW16 DM-file without checks. This may cause unspecified problems."<<endl;
	ifstream input(FileName);
	int nTemp;
	RealType fTemp;
	for(int i(0);i<m_nNum;i++)
	{
		if(input.peek()==EOF) ERROR("Length and/or shape of the YMW16 DM-file is of not compatible with the provided population");
		input>>nTemp>>fTemp>>fTemp>>fTemp>>pPulsar[i].fDM;
	}
	input.close();



}



void CPopulation::CheckVisibleRadio()
{
	for(int i(0);i<m_nNum;i++)
	{
		//cerr<<i<<" "<<pPulsar[i].GalacticPos.fL<<" "<<pPulsar[i].GalacticPos.fB<<" "<<pPulsar[i].GalacticPosDeg.fL<<" "<<pPulsar[i].GalacticPosDeg.fB<<endl;
		//cerr<<i<<" "<<m_RadioObs.CheckParkes(pPulsar[i])<<endl;	
		m_RadioObs.Print4StefansTest(&(pPulsar[i]));
	}

}

CPopulation::CPopulation()
{
	m_pCfgPtr = CConfiguration::GetInstance();
	m_pRanPtr = m_pCfgPtr->pRan;

	//m_nNum = 100;
	//m_nNum = MAXIMUM(m_pCfgPtr->Sim.nMaxAvailablePulsarNum,m_pCfgPtr->Sim.nMaxPulsarNum);
	m_nNum = m_pCfgPtr->Sim.nMaxPulsarNum;
	m_fMaxAge = m_pCfgPtr->Sim.fMaxPulsarAge;
	m_fMaxTimeStep = m_pCfgPtr->Sim.fMaxTimeStep;
//	m_pProj = new CProjection(2048,1024);

	pPulsar = new CPulsar[m_nNum];

	m_fSunAngle=0.;

	bHDFIO = false;
}

CPopulation::~CPopulation()
{
	delete [] pPulsar;
	if(bHDFIO) delete pHDFBuf;
}

RealType CPopulation::m_PulsarStartingAgeDistribution() const
{
	return m_pRanPtr->Flat() * m_fMaxAge;
}

RealType CPopulation::m_PulsarStartingAgeDistribution(const RealType fMinAge, const RealType fMaxAge) const
{
    return m_pRanPtr->Flat() * (fMaxAge-fMinAge) + fMinAge;
}


void CPopulation::Init1Per100(const bool bSmearAge,const RealType fSmearYearsSigma)
{
	const double fTimeSeparation(1e-4);//every hunderd years
	m_fMaxAge = m_nNum*fTimeSeparation;
	m_pCfgPtr->Sim.fMaxPulsarAge = m_fMaxAge;
	for(int i(0);i<m_nNum;i++)
	{
		if(bSmearAge)
		{
			pPulsar[i].fRelativeAge = - (i*fTimeSeparation + 1e-6*m_pRanPtr->Gauss(fSmearYearsSigma));
		}
		else pPulsar[i].fRelativeAge = - i*fTimeSeparation;
		pPulsar[i].Init(pPulsar[i].fRelativeAge );	
	}
}


void CPopulation::Init()
{
	for(int i(0);i<m_nNum;i++)
	{
		//pPulsar[i].fRelativeAge = -m_fMaxTimeStep*10;
		pPulsar[i].fAge=0.;
		pPulsar[i].fRelativeAge = -m_PulsarStartingAgeDistribution();
		pPulsar[i].Init(pPulsar[i].fRelativeAge );	
	}
}

void CPopulation::Init(const RealType fMinAge, const RealType fMaxAge)
{
	for(int i(0);i<m_nNum;i++)
	{
		//pPulsar[i].fRelativeAge = -m_fMaxTimeStep*10;
		pPulsar[i].fAge=0.;
		pPulsar[i].fRelativeAge = -m_PulsarStartingAgeDistribution(fMinAge,fMaxAge);
		pPulsar[i].Init(pPulsar[i].fRelativeAge );	
	}
}


void CPopulation::Zero()
{
	if(m_pCfgPtr->Sim.bSameSeedInMCMCModels)
	{
		m_pCfgPtr->Reseed( m_pCfgPtr->Model.nSeed);

	}


	if(m_pCfgPtr->Sim.bRecycleDynamics && m_pCfgPtr->Sim.bDynamicsComputed)
	{
		for(int i(0);i<m_nNum;i++)
		{
			m_pCfgPtr->Sim.bDynamics=false;
			pPulsar[i].fRelativeAge = -pPulsar[i].fAge;
			pPulsar[i].fAge=0.;
			pPulsar[i].Init(pPulsar[i].fRelativeAge );
		}	

	}
	else
	{
		Init();
	}


}


void CPopulation::Evolve()
{
	bool bCompute;
	RealType fEvoTime(0.), fLocalTimeStep(0.);

	int nSteps(0);

	if(m_pCfgPtr->Sim.bVerbose)
	{
		cout<<" Dynamics"<<endl;
		cout<<(m_pCfgPtr->Sim.bDynamicsComputed  && m_pCfgPtr->Sim.bRecycleDynamics)<<" "<<m_nNum<<endl;
	}

	for(int i(0);i<m_nNum;i++)
	{
		if(m_pCfgPtr->Sim.bVerbose && i%1000 == 0)
		{
			cout<<"  (inloop) pulsars done: "<<i<<endl;
		}
	
		//If you recycle spatial population and i-pulsar is outside the Galaxy -> do not compute
//		if((m_pCfgPtr->Sim.bDynamicsComputed  && m_pCfgPtr->Sim.bRecycleDynamics) || !pPulsar[i].bInsideGalaxy) bCompute = false;
		if( m_pCfgPtr->Sim.bDynamicsComputed  && m_pCfgPtr->Sim.bRecycleDynamics ) bCompute = false;
		else bCompute = true;

		fEvoTime= - pPulsar[i].fRelativeAge;

		if(fEvoTime > 0 && bCompute)
		{
			//obliczanie kroku czasowego niewiekszego niz maxtimestep	

			nSteps=static_cast<int>(fEvoTime/m_fMaxTimeStep)+1;
			fLocalTimeStep=fEvoTime/nSteps;

			for(int j(0);j<nSteps;j++) pPulsar[i].EvolveDynamics(fLocalTimeStep);

		}

		//musi byc przed EvolvePhysics!
		//TODO Dlaczego???

	
		pPulsar[i].EvolvePhysics(fEvoTime);

		//Note the times
		pPulsar[i].fRelativeAge = 0.;
		pPulsar[i].fAge = fEvoTime ;
		pPulsar[i].Project();
	}
	
	if(m_pCfgPtr->Sim.bVerbose) cout<<"  Pulsars done: "<<m_nNum<<endl;
	m_pCfgPtr->Sim.bDynamicsComputed=true;	


}

//do dynamicznych krokow
#if 0
void CPopulation::Evolve(const RealType fMaxTimeStep)
{
	for(int i(0);i<m_nNum;i++)
	{
		pPulsar[i].Evolve(fMaxTimeStep);
		pPulsar[i].fRelativeAge += m_fMaxTimeStep;
	}
}
#endif

void CPopulation::Project(const char *FileName)
{
	Galactic Projected;
	CProjection pProj(512,256);
	pProj.Zero();
	for(int i(0);i<m_nNum;i++)
	{
		pProj.AddPoint( (pPulsar[i].GalacticPos.fL>M_PI ?pPulsar[i].GalacticPos.fL-2*M_PI: pPulsar[i].GalacticPos.fL) ,pPulsar[i].GalacticPos.fB);
	}

	pProj.SphereDensify();
	pProj.Norm();
	pProj.Write(FileName);
}






void CPopulation::Position(const char *FileName)
{
	Cartesian Pos;
        CAreaDensity AD(-25.,-25.,25.,25.,256,256);
	for(int i(0);i<m_nNum;i++)
        {
		Pos=pPulsar[i].GetPosition();
                AD.AddPoint(Pos.fX, Pos.fY);
        }
        AD.Densify();
        AD.Norm();
        AD.Write(FileName);


}

void CPopulation::WriteBinary(const char *Filename) const
{
	ofstream out(Filename,std::ofstream::binary);
	const size_t nSize = sizeof(SBasicStar);

	SBasicStar Buffer;
	for(int i(0);i<m_nNum;i++) //if(!pPulsar[i].bIsLost) 
	{
		if(m_pCfgPtr->Sim.bVerbose && i%1000 == 0)
		{
			cout<<"  (inloop) pulsars written: "<<i<<endl;
		}
		Buffer=pPulsar[i];
		//Buffer.Data=dynamic_cast<SBasicStar&>(pPulsar[i]);
		out.write(reinterpret_cast<char *>(&Buffer),nSize);
		//out.write(static_cast<char*>(static_cast<void*>(pPulsar+i)),nSize);//aby nie uzywac reinterpret cast - stackoverflow twierdzi, ze to stabilniejsze
	}
	out.close();
}
void CPopulation::ReadBinary(const char *Filename) 
{
	ifstream in(Filename,std::ofstream::binary);
	const size_t nSize = sizeof(SBasicStar);


	//get the size of the file
	in.seekg( 0, std::ios::end );
	size_t nFileSize =  in.tellg();
	in.seekg (0, in.beg);

	int nNum=nFileSize/nSize;
	if(nFileSize%nSize > 0) ERROR("Struct size mismatch. Population's file reading not possible.");
//	nNum= MINIMUM(nNum,m_pCfgPtr->Sim.nMaxAvailablePulsarNum);
	if(nNum>m_pCfgPtr->Sim.nMaxPulsarNum) WARNING("Input larger than maximum allowed population. Reading only up to max.");
	nNum= MINIMUM(nNum,m_pCfgPtr->Sim.nMaxPulsarNum);
	m_pCfgPtr->Sim.nMaxPulsarNum=m_nNum=nNum;

//	if(pPulsar!=NULL) delete [] pPulsar;
//	pPulsar = new CPulsar[m_nNum];

	SBasicStar Buffer;
	if(	m_pCfgPtr->Sim.bVerbose ) cout<<"Num: "<<nNum<<endl;
	//in.read(static_cast<char*>(static_cast<void*>(pPulsar)),nSize*nNum);
	for(int i(0);i<nNum;i++)
	{
		in.read(reinterpret_cast<char *>(&Buffer),nSize);
		pPulsar[i]=Buffer;
	}
	m_pCfgPtr->Sim.bDynamics=false;//FIXME!!!!!!!!! szybki hack 18.05.15
	m_pCfgPtr->Sim.bDynamicsComputed=true;
	in.close();


}

void CPopulation::ReadBinaryRelevantGeometry(const char *Filename) 
{
	ifstream in(Filename,std::ofstream::binary);
	const size_t nSize = sizeof(SBasicStar);


	//get the size of the file
	in.seekg( 0, std::ios::end );
	size_t nFileSize =  in.tellg();
	in.seekg (0, in.beg);

	int nNum=nFileSize/nSize;
	if(nFileSize%nSize > 0) ERROR("Struct size mismatch. Population's file reading not possible.");
	//nNum= MINIMUM(nNum,m_pCfgPtr->Sim.nMaxPulsarNum);


	CSurvey Srv;

	int nCorrect(0);
	SBasicStar Buffer;
	for(int i(0);(i<nNum && nCorrect<m_pCfgPtr->Sim.nMaxPulsarNum);i++)
	{
		in.read(reinterpret_cast<char *>(&Buffer),nSize);
		if( Srv.IsCoveredInSurvey(Buffer.GalacticPosDeg.fL,Buffer.GalacticPosDeg.fB, ParkesMultiBeam) && Buffer.GalacticPosDeg.fDist <= m_pCfgPtr->Sim.fMaxDist  )
		{
			pPulsar[nCorrect]=Buffer;
			nCorrect++;
		}
	}
	cout<<"nNum: "<<nNum<<" nCorrect: "<<nCorrect<<endl;

	m_pCfgPtr->Sim.nMaxPulsarNum=m_nNum=nCorrect;

	m_pCfgPtr->Sim.bDynamics=false;//FIXME!!!!!!!!! szybki hack 18.05.15
	m_pCfgPtr->Sim.bDynamicsComputed=true;
	in.close();

}



void CPopulation::Write(const char *Filename) const
{
	ofstream out(Filename);

	//CAreaDensity logPPdot(m_pCfgPtr->Plot.fZeroLogP,m_pCfgPtr->Plot.fZeroLogPdot,m_pCfgPtr->Plot.fMaxLogP,m_pCfgPtr->Plot.fMaxLogPdot, m_pCfgPtr->Plot.nResLogP, m_pCfgPtr->Plot.nResLogPdot);
	//CAreaDensity logPPdotParkes(m_pCfgPtr->Plot.fZeroLogP,m_pCfgPtr->Plot.fZeroLogPdot,m_pCfgPtr->Plot.fMaxLogP,m_pCfgPtr->Plot.fMaxLogPdot, m_pCfgPtr->Plot.nResLogP, m_pCfgPtr->Plot.nResLogPdot);


	//results' description:	
	out<<"# Age[Myr] X Y Z [kpc] VX VY VZ [kpc/Myr] L B [deg] Dist[kpc]";

	if(m_pCfgPtr->Sim.bPhysics)
	{
		out<<" Period[s] PeriodDot[s/s] BField[Gauss]";
	}

	if(m_pCfgPtr->Sim.bObservations)
	{
		out<<" DM[pc/(cm^3)] L400[mJy*(kpc^2)] VisibleInParkes[bool]";
	}

	out<<endl;


	out<<setprecision(14);
	for(int i(0);i<m_nNum;i++)
	{
		if(! pPulsar[i].bIsLost)
		{	
			out<<scientific;
			out<<pPulsar[i].fAge<<" ";
			out<<pPulsar[i].Position.fX<<" "<<pPulsar[i].Position.fY<<" "<<pPulsar[i].Position.fZ<<" ";
			out<<pPulsar[i].Velocity.fX<<" "<<pPulsar[i].Velocity.fY<<" "<<pPulsar[i].Velocity.fZ<<" ";
			out<<pPulsar[i].GalacticPosDeg.fL<<" "<<pPulsar[i].GalacticPosDeg.fB<<" "<<pPulsar[i].GalacticPosDeg.fDist;

			if(m_pCfgPtr->Sim.bPhysics)
			{
				out<<" "<<pPulsar[i].fPeriod<<" "<<pPulsar[i].fPeriodDot<<" "<<pPulsar[i].fBField;

				if(m_pCfgPtr->Sim.bObservations)
				{
					out<<" "<<pPulsar[i].fDM<<" "<<pPulsar[i].fL400<<" "<<m_RadioObs.CheckParkes(&(pPulsar[i]));
				}
			}
			out<<endl;
		}
	}
	out.close();

}


void CTestPopulation::WriteReadTest()
{
	
	Init();
TTimestamp T0=GetTimestamp();
	Evolve();
TTimestamp T1=GetTimestamp();
	WriteBinary("Out.bin");
	Write("A.dat");
	Zero();
TTimestamp T2=GetTimestamp();

	Evolve();
TTimestamp T3=GetTimestamp();

	WriteBinary("Out2.bin");
	Write("B.dat");

	cout<<(T1-T0)/1e6<<" "<<(T3-T2)/1e6<<endl;
}



#if 0
Powienien robic test
void CTestPopulation::EnergyTest()
{

	

	Init();

	ofstream out("CTestPopulation.Energy.dat");


	RealType fMaxTime(250.),fDTime(0.0001);
	
	RealType fEKin(0.), fEPot(0.), fETot(0.);
	RealType fMass(1.), fVTot(0.);
	Cartesian L;
RealType fETotBeg, fLZBeg;
	out<<scientific<<setprecision(19);
	for(int i(0);i<=25000;i++)
	{	
		fVTot=Length(P1.m_Velocity);
		fEKin=fMass*fVTot*fVTot*0.5;
		fEPot=fMass*P1.GalacticPotential(P1.m_Position);
		fETot=fEKin - fEPot;
		L=1.*CrossProduct(P1.m_Position,P1.m_Velocity);//masa 1.

		(P1.m_Position);
		//P1.m_EvolveRKi3Adaptive(0.1);
		out<<i*0.001<<" "<<P1.m_Position.fX<<" "<<P1.m_Position.fY<<" "<<P1.m_Position.fZ<<" "<<fEKin<<" "<<-fEPot<<" "<<fETot<<" "<<L.fZ<<endl;
		
		if(i==0)
		{
			fETotBeg=fETot;
			fLZBeg=L.fZ;
		}


		P1.Evolve(0.001);
	}
	cout<<scientific<<setprecision(19)<<fETot-fETotBeg<<" "<<L.fZ-fLZBeg<<endl;

	out.close();
}
#endif



void CPopulation::WriteEvolutionTracksHDF5(const char *FileName, const int nStepsPerTrack, const int *pNumery, const int nNumNum)
{
	CPulsar *pTrack = new CPulsar[nStepsPerTrack+2];
	CPulsar Current;

	double *pfBuf = new double[m_nNum];
	
	// Create HDF5 file and dataset
	H5::H5File HDF5File(FileName, H5F_ACC_TRUNC);
	const int nRank(1);
	hsize_t pnDims[nRank];
	//size of single track
	pnDims[0]=static_cast<hsize_t>(nStepsPerTrack+2);
    	H5::DataSpace Space(nRank, pnDims);

	H5::DataSet Set;
	H5::Group Group;

	//set to recompute only physics
	m_pCfgPtr->Sim.bRecycleDynamics=true;
	m_pCfgPtr->Sim.bDynamics=false;
	m_pCfgPtr->Sim.bPhysics=true;
	m_pCfgPtr->Sim.bObservations=true;

	RealType fDT,fMaxAge;//in [Myr]
	stringstream Translator;
	//loop for the whole pop
	for(int k(0);k<nNumNum;k++)
	//for(int k(0);k<1000;k++)
	{

		//Init Single Star
		Current=pPulsar[pNumery[k]];
		if(m_pCfgPtr->Sim.bVerbose) 
		{
			Current.PrintDescriptiveStarInfo(cout);
		}
		//reset age to initial values

		Current.ResetToInit(m_pCfgPtr->Sim.bPhysics,m_pCfgPtr->Sim.bDynamics);

		//Current.fMinimalBField=pow(10., 8. );//FIXME!!!


		//save the init state
		pTrack[0]=Current;

		//Compute dt
		fDT=-Current.fRelativeAge /(nStepsPerTrack+1);		
		if(m_pCfgPtr->Sim.bVerbose) 
		{
			cout<<Current.fPeriod<<" "<<Current.fPeriodDot<<" "<<Current.fBField<<" "<<fDT<<" "<<Current.fAge<<endl;
		}
		//Evolve it
		for(int j(1);j<=nStepsPerTrack;j++)
		{
			Current.EvolvePhysics(fDT);
			Current.fAge+=fDT;
			Current.fRelativeAge+=fDT;
			pTrack[j]=Current;
		}

		//Evolve in one step
		Current=pTrack[0];//the zero step
		//Current.fMinimalBField=pow(10., 8. );//FIXME!!!
		Current.CPulsarPhysics::Evolve(fDT*(nStepsPerTrack+1),true);
		pTrack[nStepsPerTrack+1]=Current;

	  
		Translator.str("");
		Translator.clear();
		Translator<<"/"<<k;
		//cout<<Translator.str().c_str()<<endl;

		Group = HDF5File.createGroup( Translator.str().c_str() );
		//Write track to the HDF
		
		Translator.str("");
		Translator.clear();
		Translator<<"/"<<k<<"/fPeriod";
		//cout<<Translator.str().c_str()<<endl;

		Set = HDF5File.createDataSet(Translator.str().c_str(), H5::PredType::NATIVE_DOUBLE, Space);
		for(int i(nStepsPerTrack+2);i--;) pfBuf[i]=pTrack[i].fPeriod;
		Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

		Translator.str("");
		Translator.clear();
		Translator<<"/"<<k<<"/fPeriodDot";
		//cout<<Translator.str().c_str()<<endl;

		Set = HDF5File.createDataSet(Translator.str().c_str(), H5::PredType::NATIVE_DOUBLE, Space);
		for(int i(nStepsPerTrack+2);i--;) pfBuf[i]=pTrack[i].fPeriodDot;
		Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

		Translator.str("");
		Translator.clear();
		Translator<<"/"<<k<<"/fBField";
		//cout<<Translator.str().c_str()<<endl;

		Set = HDF5File.createDataSet(Translator.str().c_str(), H5::PredType::NATIVE_DOUBLE, Space);
		for(int i(nStepsPerTrack+2);i--;) pfBuf[i]=pTrack[i].fBField;
		Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);
		if(m_pCfgPtr->Sim.bVerbose && k%100 == 0)
		{
			cout<<"  (inloop) pulsars' tracks done: "<<k<<endl;
		}
	}


	delete [] pTrack;
	delete [] pfBuf;

}






CHDFPopulationBuffer::CHDFPopulationBuffer(const char *FileName)
{
	m_pCfgPtr = CConfiguration::GetInstance();

	pFile = new H5::H5File( FileName, H5F_ACC_TRUNC );

	pDims = new hsize_t[cnDims];
	pDimsMax = new hsize_t[cnDims];
	pOffset = new hsize_t[cnDims];


	pDims[0]=static_cast<hsize_t>(m_pCfgPtr->Sim.nMaxPulsarNum);
	pDimsMax[0]=static_cast<hsize_t>(m_pCfgPtr->Sim.nMaxPulsarNum * m_pCfgPtr->Sim.nNumberOfIteration);

	pbBuf = new short[pDims[0]];
	pfBuf = new double[pDims[0]];	

	CreateDatasets();
}

CHDFPopulationBuffer::~CHDFPopulationBuffer()
{
	delete [] pDims;
	delete [] pDimsMax;
	delete [] pOffset;
	delete [] pbBuf;
	delete [] pfBuf;
}

void CHDFPopulationBuffer::CreateDatasets()
{
	if(m_pCfgPtr->Sim.bVerbose) cout<<"Creating datasets"<<endl;

	H5::DataSpace DataSpace(cnDims, pDimsMax, pDimsMax);

  	H5::DSetCreatPropList Plist;
	Plist.setChunk(cnDims,pDims);

	Set.bBeaming = pFile->createDataSet("bBeaming", H5::PredType::NATIVE_SHORT, DataSpace,Plist );
	Set.bInsideGalaxy = pFile->createDataSet("bInsideGalaxy", H5::PredType::NATIVE_SHORT, DataSpace, Plist);
	Set.bIsLost = pFile->createDataSet("bIsLost", H5::PredType::NATIVE_SHORT, DataSpace, Plist);
	Set.fProbability = pFile->createDataSet("fProbability", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fAge = pFile->createDataSet("fAge", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fRelativeAge = pFile->createDataSet("fRelativeAge", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fMass = pFile->createDataSet("fMass", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fPeriod = pFile->createDataSet("fPeriod", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fPeriodDot = pFile->createDataSet("fPeriodDot", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fBField = pFile->createDataSet("fBField", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialPeriod = pFile->createDataSet("fInitialPeriod", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialPeriodDot = pFile->createDataSet("fInitialPeriodDot", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialBField = pFile->createDataSet("fInitialBField", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fSpectralIndex = pFile->createDataSet("fSpectralIndex", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fL400 = pFile->createDataSet("fL400", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fS400 = pFile->createDataSet("fS400", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fS1400 = pFile->createDataSet("fS1400", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fDM = pFile->createDataSet("fDM", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fPositionX = pFile->createDataSet("fPositionX", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fPositionY = pFile->createDataSet("fPositionY", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fPositionZ = pFile->createDataSet("fPositionZ", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fVelocityX = pFile->createDataSet("fVelocityX", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fVelocityY = pFile->createDataSet("fVelocityY", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fVelocityZ = pFile->createDataSet("fVelocityZ", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialPositionX = pFile->createDataSet("fInitialPositionX", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialPositionY = pFile->createDataSet("fInitialPositionY", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialPositionZ = pFile->createDataSet("fInitialPositionZ", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialKickX = pFile->createDataSet("fInitialKickX", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialKickY = pFile->createDataSet("fInitialKickY", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fInitialKickZ = pFile->createDataSet("fInitialKickZ", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fGalacticPosB = pFile->createDataSet("fGalacticPosB", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fGalacticPosL = pFile->createDataSet("fGalacticPosL", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fGalacticPosDist = pFile->createDataSet("fGalacticPosDist", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fGalacticPosDegB = pFile->createDataSet("fGalacticPosDegB", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fGalacticPosDegL = pFile->createDataSet("fGalacticPosDegL", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.fBeamingFraction = pFile->createDataSet("fBeamingFraction", H5::PredType::NATIVE_DOUBLE, DataSpace, Plist);
	Set.bVisibleInParkes = pFile->createDataSet("bVisibleInParkes", H5::PredType::NATIVE_SHORT, DataSpace, Plist);
	Set.bVisibleInSKA = pFile->createDataSet("bVisibleInSKA", H5::PredType::NATIVE_SHORT, DataSpace, Plist);
	Set.bVisibleInSKALC = pFile->createDataSet("bVisibleInSKALC", H5::PredType::NATIVE_SHORT, DataSpace, Plist);
	if(m_pCfgPtr->Sim.bVerbose) cout<<"Datasets created"<<endl;

}


void CHDFPopulationBuffer::WriteChunk(const CPulsar *pPulsar, const int nIteration)
{
	H5::DataSpace SlabSpace;
	pOffset[0]=nIteration*pDims[0];

	//Turned out that the example was wrong the last argument can't be larger then available space + offset -> generates exeption	
	H5::DataSpace DataSpace(cnDims, pDims, pDims);

	if(m_pCfgPtr->Sim.bVerbose)
	{
		cout<<"Writing a chunk of data."<<endl;
		cout<<"nIteration: "<<nIteration<<endl;
		cout<<"offset: "<<pOffset[0]<<endl;
	}

	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pbBuf[i]=pPulsar[i].bBeaming;
	SlabSpace = Set.bBeaming.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.bBeaming.write( pbBuf, H5::PredType::NATIVE_SHORT, DataSpace, SlabSpace );

	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pbBuf[i]=pPulsar[i].bInsideGalaxy;
	SlabSpace = Set.bInsideGalaxy.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.bInsideGalaxy.write( pbBuf, H5::PredType::NATIVE_SHORT, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pbBuf[i]=pPulsar[i].bIsLost;
	SlabSpace = Set.bIsLost.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.bIsLost.write( pbBuf, H5::PredType::NATIVE_SHORT, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fProbability;
	SlabSpace = Set.fProbability.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fProbability.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fAge;
	SlabSpace = Set.fAge.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fAge.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fRelativeAge;
	SlabSpace = Set.fRelativeAge.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fRelativeAge.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fMass;
	SlabSpace = Set.fMass.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fMass.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fPeriod;
	SlabSpace = Set.fPeriod.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fPeriod.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fPeriodDot;
	SlabSpace = Set.fPeriodDot.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fPeriodDot.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fBField;
	SlabSpace = Set.fBField.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fBField.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fInitialPeriod;
	SlabSpace = Set.fInitialPeriod.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialPeriod.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fInitialPeriodDot;
	SlabSpace = Set.fInitialPeriodDot.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialPeriodDot.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fInitialBField;
	SlabSpace = Set.fInitialBField.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialBField.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fSpectralIndex;
	SlabSpace = Set.fSpectralIndex.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fSpectralIndex.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fL400;
	SlabSpace = Set.fL400.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fL400.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fS400;
	SlabSpace = Set.fS400.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fS400.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fS1400;
	SlabSpace = Set.fS1400.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fS1400.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].fDM;
	SlabSpace = Set.fDM.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fDM.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );

	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].Position.fX;
	SlabSpace = Set.fPositionX.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fPositionX.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].Position.fY;
	SlabSpace = Set.fPositionY.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fPositionY.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].Position.fZ;
	SlabSpace = Set.fPositionZ.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fPositionZ.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].Velocity.fX;
	SlabSpace = Set.fVelocityX.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fVelocityX.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].Velocity.fY;
	SlabSpace = Set.fVelocityY.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fVelocityY.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].Velocity.fZ;
	SlabSpace = Set.fVelocityZ.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fVelocityZ.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].InitialPosition.fX;
	SlabSpace = Set.fInitialPositionX.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialPositionX.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].InitialPosition.fY;
	SlabSpace = Set.fInitialPositionY.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialPositionY.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].InitialPosition.fZ;
	SlabSpace = Set.fInitialPositionZ.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialPositionZ.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].InitialKick.fX;
	SlabSpace = Set.fInitialKickX.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialKickX.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].InitialKick.fY;
	SlabSpace = Set.fInitialKickY.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialKickY.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );

	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].InitialKick.fZ;
	SlabSpace = Set.fInitialKickZ.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fInitialKickZ.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].GalacticPos.fB;
	SlabSpace = Set.fGalacticPosB.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fGalacticPosB.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].GalacticPos.fL;
	SlabSpace = Set.fGalacticPosL.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fGalacticPosL.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].GalacticPos.fDist;
	SlabSpace = Set.fGalacticPosDist.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fGalacticPosDist.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].GalacticPosDeg.fB;
	SlabSpace = Set.fGalacticPosDegB.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fGalacticPosDegB.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );
	
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) pfBuf[i]=pPulsar[i].GalacticPosDeg.fL;
	SlabSpace = Set.fGalacticPosDegL.getSpace();
	SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
	Set.fGalacticPosDegL.write( pfBuf, H5::PredType::NATIVE_DOUBLE, DataSpace, SlabSpace );


	for(int s(0);s<m_pCfgPtr->SurveySelection.nNum;s++)
	{
		if(m_pCfgPtr->SurveySelection.pID[s] == ParkesMultiBeam )
		{
			for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) 
			{
				if( !pPulsar[i].bIsLost && (pPulsar[i].fBField >= 1e10) ) pbBuf[i]=m_RadioObs.IsVisible(&(pPulsar[i]), ParkesMultiBeam);
				else pbBuf[i]=false;
			}
			SlabSpace = Set.bVisibleInParkes.getSpace();
			SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
			Set.bVisibleInParkes.write( pbBuf, H5::PredType::NATIVE_SHORT, DataSpace, SlabSpace );
		}
		if(m_pCfgPtr->SurveySelection.pID[s] == SKA1Mid )
		{
			for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) 
			{
				if( !pPulsar[i].bIsLost && (pPulsar[i].fBField >= 1e10) ) pbBuf[i]=m_RadioObs.IsVisible(&(pPulsar[i]), SKA1Mid);
				else pbBuf[i]=false;
			}
			SlabSpace = Set.bVisibleInSKA.getSpace();
			SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
			Set.bVisibleInSKA.write( pbBuf, H5::PredType::NATIVE_SHORT, DataSpace, SlabSpace );
		}
		if(m_pCfgPtr->SurveySelection.pID[s] == SKA1MidLC )
		{
			for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) 
			{
				if( !pPulsar[i].bIsLost && (pPulsar[i].fBField >= 1e10) ) pbBuf[i]=m_RadioObs.IsVisible(&(pPulsar[i]), SKA1MidLC);
				else pbBuf[i]=false;
			}
			SlabSpace = Set.bVisibleInSKALC.getSpace();
			SlabSpace.selectHyperslab( H5S_SELECT_SET, pDims, pOffset );
			Set.bVisibleInSKALC.write( pbBuf, H5::PredType::NATIVE_SHORT, DataSpace, SlabSpace );
		}



	}	

}

void CPopulation::WriteHDF5(const char *FileName, const int nIteration)
{
	if(!bHDFIO) 
	{
		pHDFBuf = new CHDFPopulationBuffer(FileName);
		bHDFIO=true;
	}


	pHDFBuf->WriteChunk(pPulsar,nIteration);

}

void CPopulation::WriteLost(const char* FileName)
{
	ofstream Out(FileName);
	for(int i(m_pCfgPtr->Sim.nMaxPulsarNum);i--;) if( pPulsar[i].bIsLost )
	{
		Out<<i<<" "<<pPulsar[i].ETot()<<" "<<pPulsar[i].fAge<<" ";
		Out<<pPulsar[i].Position.fX<<" "<<pPulsar[i].Position.fY<<" "<<pPulsar[i].Position.fZ<<" ";
		Out<<pPulsar[i].Velocity.fX<<" "<<pPulsar[i].Velocity.fY<<" "<<pPulsar[i].Velocity.fZ<<" ";
		Out<<endl;
	}
	Out.close();

}



#if 0
void CPopulation::WriteHDF5(const char *FileName) const
{

//NATIVE_DOUBLE - to jest double
//
//NATIVE_SHORT unsigned short int
//
    // Data to write

    // the length of the data
//    int length = sizeof(person_list) / sizeof(PersonalInformation);
    
// the array of each length of multidimentional data.
//    dim[0] = sizeof(person_list) / sizeof(PersonalInformation);
    //dim[0] = nNum;

    // the length of dim
//sizeof(dim) / sizeof(hsize_t);///FIXME czym to jest?!
	short *pbBuf = new short[m_nNum];
	double *pfBuf = new double[m_nNum];

   // Create HDF5 file and dataset
	H5::H5File HDF5File(FileName, H5F_ACC_TRUNC);
	const int nRank(1);
	hsize_t pnDims[nRank];
	pnDims[0]=static_cast<hsize_t>(m_nNum);
    	H5::DataSpace Space(nRank, pnDims);

	H5::DataSet Set;

	Set = HDF5File.createDataSet("bBeaming", H5::PredType::NATIVE_SHORT, Space);
	for(int i(m_nNum);i--;) pbBuf[i]=pPulsar[i].bBeaming; //Rearenge the data
	Set.write(pbBuf,H5::PredType::NATIVE_SHORT);

	Set = HDF5File.createDataSet("bInsideGalaxy", H5::PredType::NATIVE_SHORT, Space);
	for(int i(m_nNum);i--;) pbBuf[i]=pPulsar[i].bInsideGalaxy; //Rearenge the data
	Set.write(pbBuf,H5::PredType::NATIVE_SHORT);

	Set = HDF5File.createDataSet("bIsLost", H5::PredType::NATIVE_SHORT, Space);
	for(int i(m_nNum);i--;) pbBuf[i]=pPulsar[i].bIsLost; //Rearenge the data
	Set.write(pbBuf,H5::PredType::NATIVE_SHORT);



	Set = HDF5File.createDataSet("fProbability", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fProbability; //Rearenge the data
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);


	Set = HDF5File.createDataSet("fAge", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fAge;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fRelativeAge", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fRelativeAge;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fMass", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fMass;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fPeriod", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fPeriod;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fPeriodDot", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fPeriodDot;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fBField", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fBField;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fInitialPeriod", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fInitialPeriod;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fInitialPeriodDot", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fInitialPeriodDot;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fInitialBField", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fInitialBField;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fSpectralIndex", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fSpectralIndex;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fL400", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fL400;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fS400", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fS400;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fS1400", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fS1400;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fDM", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].fDM;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);



	Set = HDF5File.createDataSet("fPositionX", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].Position.fX;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fPositionY", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].Position.fY;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fPositionZ", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].Position.fZ;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fVelocityX", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].Velocity.fX;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fVelocityY", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].Velocity.fY;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fVelocityZ", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].Velocity.fZ;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);


	Set = HDF5File.createDataSet("fInitialPositionX", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].InitialPosition.fX;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fInitialPositionY", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].InitialPosition.fY;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fInitialPositionZ", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].InitialPosition.fZ;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);


	Set = HDF5File.createDataSet("fInitialKickX", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].InitialKick.fX;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fInitialKickY", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].InitialKick.fY;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fInitialKickZ", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].InitialKick.fZ;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fGalacticPosB", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].GalacticPos.fB;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fGalacticPosL", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].GalacticPos.fL;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fGalacticPosDist", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].GalacticPos.fDist;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);



	Set = HDF5File.createDataSet("fGalacticPosDegB", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].GalacticPosDeg.fB;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fGalacticPosDegL", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].GalacticPosDeg.fL;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	//this one is redundant
//	Set = HDF5File.createDataSet("fGalacticPosDegDist", H5::PredType::NATIVE_DOUBLE, Space);
//	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].GalacticPosDeg.fDist;
//	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("bVisibleInParkes", H5::PredType::NATIVE_SHORT, Space);
	for(int i(m_nNum);i--;) pbBuf[i]=m_RadioObs.CheckParkes(pPulsar[i]) && !pPulsar[i].bIsLost;
	Set.write(pbBuf,H5::PredType::NATIVE_SHORT);

	Set = HDF5File.createDataSet("fBeamingFraction", H5::PredType::NATIVE_DOUBLE, Space);
	for(int i(m_nNum);i--;) pfBuf[i]=pPulsar[i].BeamingFactor() ;
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);


    

	delete [] pbBuf;
	delete [] pfBuf;

}

#endif


