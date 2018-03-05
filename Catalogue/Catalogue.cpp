#include<fstream>
#include<iostream>
#include<algorithm>


#include<string>
#include<sstream> 
#include<cmath>
#include<limits>
#include<bitset>

#include<Catalogue.hpp>
#include<NANcinHandler.hpp>

#include<Debug.hpp>

using namespace std;


#include<ModelParameters.hpp>
#include<Pulsar.hpp>

SSurveys SurveyID2Word(ESurveyID ID)
{
	SSurveys Temp;
	Temp.Word=0;

	//FIXME zamienic na switch-case
	if(ID == ParkesMultiBeam )
	{
		Temp.Mask.bParkesMultiBeam=true;
	}
	if(ID == Parkes70 ) 
	{
		Temp.Mask.bParkes70=true;
	}
	if(ID == SwinIntermLat)
	{
		Temp.Mask.bSwinIntermLat=true;
	}
	if(ID == SwinExtended )
	{
		Temp.Mask.bSwinExtended=true;
	}


	return Temp;

}

CCatalogue::CCatalogue(const char * FileName)
{

	m_pCfgPtr=CConfiguration::GetInstance();
	ReadTxt(FileName);
	//m_pCfgPtr->pATNF=this;
}


CCatalogue::CCatalogue()
{

	m_pCfgPtr=CConfiguration::GetInstance();
	ReadTxt(m_pCfgPtr->Sim.sCatalogueFile);
    pbUsed = new bool[nTotalNum];
    for(int i(0);i<nTotalNum;i++) pbUsed[i]=false;
	//m_pCfgPtr->pATNF=this;
}





CCatalogue::~CCatalogue()
{
	delete [] pPulsar;
    delete [] pbUsed;
}


SSurveys CCatalogue::Check4Surveys(const string Input) const
{
	SSurveys Temp;
	Temp.Word=0;
	//FIXME zaimplemenetować resztę przeglądów o których napisał mi Stefan
	if(Input.find("pksmb") != string::npos)
	{
		Temp.Mask.bParkesMultiBeam=true;
	}
	if(Input.find("pks70") !=  string::npos)
	{
		Temp.Mask.bParkes70=true;
	}
	if(Input.find("pkssw") !=  string::npos)
	{
		Temp.Mask.bSwinIntermLat=true;
	}
	if(Input.find("pkshl") !=  string::npos)
	{
		Temp.Mask.bSwinExtended=true;
	}
//	cout<<std::bitset<8>(Temp.Word)<<" "<<Input<<endl;
	return Temp;
}


void CCatalogue::ReadTxt(const char * FileName)
{
	fstream InFile(FileName);
	if(!InFile.good())
	{
		string msg("Could not open catalogue file ");
		msg+=FileName;
		ERROR(msg);
	}
	//count number of lines
	nTotalNum = count(istreambuf_iterator<char>(InFile), istreambuf_iterator<char>(), '\n');
	InFile.seekg(0);
	//cout<<"Catalogue size: "<<nTotalNum<<endl;


	pPulsar = new SCatalogueEntry[nTotalNum];
	string Nothing, Value;

	for(int i(0);i<nTotalNum;i++)
	{	
		InFile	>> Nothing;//1 numer porzadkowy

		InFile	>> pPulsar[i].Name;//2 nazwa

		InFile	>> Nothing;//3 ref
		
		InFile	>> Nothing;//4 jeszcze raz nazwa
		InFile	>> Nothing;//5 ref

		InFile	>> Value;//6 GalL
		pPulsar[i].fDegGalL = Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bDegGalL = IsPresent(Value); 

		InFile	>> Value;//7 GalB
		pPulsar[i].fDegGalB = Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bDegGalB = IsPresent(Value); 

		InFile	>> Value;//8 P0
		pPulsar[i].fP=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bP = IsPresent(Value); 


		InFile	>> Nothing;//9 dP0
		InFile	>> Nothing;//10 ref

		InFile	>> Value;//11 P1
		pPulsar[i].fPdot=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bPDot = IsPresent(Value); 
		if(pPulsar[i].fPdot<0.)
		{
			pPulsar[i].fPdot=EmitNan();
			pPulsar[i].DataFields.Mask.bPDot = false;
		}

		InFile	>> Nothing;//12 dP1
		InFile	>> Nothing;//13 ref

		InFile	>> Value;//14 DM
		pPulsar[i].fDM=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bDM = IsPresent(Value); 

		InFile	>> Nothing;//15 dDM
		InFile	>> Nothing;//16 ref

		InFile	>> Value;//17 RM
		pPulsar[i].fRM=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bRM = IsPresent(Value); 


		InFile	>> Nothing;//18 dRM
		InFile	>> Nothing;//19 ref

		InFile	>> Nothing; //20 binary
		InFile	>> Nothing; //21 ref

		InFile	>> Value;//22 Dist
		pPulsar[i].fDist=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bDist = IsPresent(Value); 

		InFile	>> Value;//23 DistDM
		pPulsar[i].fDistDM=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bDistDM = IsPresent(Value); 
		
		InFile	>> Nothing; //24 ref

		InFile	>> Value;//25 Survey
		pPulsar[i].Surveys = Check4Surveys(Value);

		InFile	>> Value;//26 Survey mask
		pPulsar[i].nSurveys=Check4Nan(Value);

		InFile	>> Value;//27 Radio Flux @ 400Mhz [mJy]
		pPulsar[i].fRadioFlux400=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bRadioFlux400 = IsPresent(Value); 
	
		InFile	>> Nothing; //28 dS400
		InFile	>> Nothing; //29 ref

		InFile	>> Value;//30 Radio Flux @ 1400Mhz [mJy]
		pPulsar[i].fRadioFlux1400=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bRadioFlux1400 = IsPresent(Value); 

		InFile	>> Nothing; //31 dS1400
		InFile	>> Nothing; //32 ref

		InFile	>> Value;//33 Age [yr]
		pPulsar[i].fAge=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bAge = IsPresent(Value); 

		InFile	>> Value;//34 BSurf[G]
		pPulsar[i].fBSurface=Check4Nan(Value);
		pPulsar[i].DataFields.Mask.bBSurface = IsPresent(Value); 

		InFile	>> Nothing;//35 P1_I

		InFile	>> Nothing;//36 Age_I

		InFile	>> Nothing;//37 Bsurf_I

	}


	InFile.close();
}


void CCatalogue::Print()
{
	std::bitset<16> x;
	std::bitset<8> y;
	bool bUsable;

	SDataFields Need;
	Need.Mask.bDegGalL=1;
	Need.Mask.bDegGalB=1;
	Need.Mask.bP=1;
	Need.Mask.bPDot=1;
	Need.Mask.bDM=1;
	Need.Mask.bRadioFlux1400=1;

	for(int i(0);i<nTotalNum;i++)
	{
		cout	<< pPulsar[i].Name <<" "
			<< pPulsar[i].fDegGalL <<" "
			<< pPulsar[i].fDegGalB <<" "
			<< pPulsar[i].fP <<" "
			<< pPulsar[i].fPdot <<" "
			<< pPulsar[i].fDM <<" "
//			<< pPulsar[i].fRM<<" "
//			<< pPulsar[i].fDist<<" "
			<< pPulsar[i].fDistDM<<" "
			<< pPulsar[i].nSurveys<<" "
//			<< pPulsar[i].fRadioLuminosity<<" "
//			<< pPulsar[i].fRadioLuminosity1400<<" "
//			<< pPulsar[i].fAge<<" "
			<< pPulsar[i].fBSurface<<" ";
		x=pPulsar[i].DataFields.Word;//nie moge tego castem zrobic??
		cout<<x<<" ";
		y=pPulsar[i].Surveys.Word;
		cout<<y<<" ";
		bUsable = (pPulsar[i].DataFields.Word & Need.Word ) == Need.Word ;
		cout<<bUsable;
		cout<<endl;
	}

}

bool CObservations::m_AboveDeathLines(const SCatalogueEntry Pulsar) const
{
	CPulsar Dummy;
	if (Dummy.AboveDeathLines(Pulsar.fP,Pulsar.fPdot)) return true;
	return false;
}
bool CObservations::m_InsideArea(const SCatalogueEntry Pulsar, const ESurveyID ID) const
{
	CSurvey Srv;
	return Srv.IsCoveredInSurvey( Pulsar.fDegGalL, Pulsar.fDegGalB, ID );

}

CObservations::CObservations(const ESurveyID ID)
{
	Vector.Zero();
	nUsedATNFPulsars=0;
	SurveyID=ID;
	SurveyMask=SurveyID2Word(SurveyID);
	Survey.Init(SurveyID);
	CConfiguration *pCfgPtr=CConfiguration::GetInstance();

	m_pCfgPtr = CConfiguration::GetInstance();

	
	ofstream ATNFCheck;
	if(pCfgPtr->Sim.bVerbose) 
	{
		ATNFCheck.open("ATNF.Check.dat");
		ATNFCheck<<"#L B P Pdot S1400 DM"<<endl;
	}
	bool bAboveDeath, bInsideArea;
	int nDead(0);
	//Not optimal but since all CObservations are constructed only once.
	//Count the number PSR's covered by given survey

	bool *pMaska = new bool[pCfgPtr->pATNF->nTotalNum];
	//cout<<pCfgPtr->pATNF->nTotalNum<<endl;
	for(int i(0); i<pCfgPtr->pATNF->nTotalNum; i++)
	{
		pMaska[i]=false;
		//cout<<"i: "<<i<<endl;
		//cout<<std::bitset<8>(pCfgPtr->pATNF->pPulsar[i].Surveys.Word & SurveyMask.Word)<<" "<<std::bitset<8>(SurveyMask.Word) <<endl;
		if((pCfgPtr->pATNF->pPulsar[i].Surveys.Word & SurveyMask.Word) == SurveyMask.Word )
		{
			if( (pCfgPtr->pATNF->pPulsar[i].DataFields.Word & pCfgPtr->Sim.DataFields.Word) == pCfgPtr->Sim.DataFields.Word)
			{
				if(pCfgPtr->Sim.bATNFDeathLinesCut && !m_AboveDeathLines(pCfgPtr->pATNF->pPulsar[i]) )
				{
					bAboveDeath=false;
					//nDead++;
				}
				else bAboveDeath=true;

				if(pCfgPtr->Sim.bCutATNFOutsideArea && !m_InsideArea(pCfgPtr->pATNF->pPulsar[i],ID))
				{
					bInsideArea = false;
				}
				else bInsideArea = true;

				if (bAboveDeath && bInsideArea)
				{
					nUsedATNFPulsars++;
					m_AddToMap(  ComputeSpacePointCoordiantesRadio( pCfgPtr->pATNF->pPulsar[i] )  );
					m_AddToMap(  ComputeSpacePointCoordiantesSpace( pCfgPtr->pATNF->pPulsar[i] )  );
					m_AddToMap(  ComputeSpacePointCoordiantes4D( pCfgPtr->pATNF->pPulsar[i] ) );

					m_AddToCM(pCfgPtr->pATNF->pPulsar[i] );

					Vector.AddPoint(pCfgPtr->pATNF->pPulsar[i] );

                    pCfgPtr->pATNF->pbUsed[i]=true;
					if(pCfgPtr->Sim.bVerbose) 
					{
						ATNFCheck<<pCfgPtr->pATNF->pPulsar[i].fDegGalL<<" ";
						ATNFCheck<<pCfgPtr->pATNF->pPulsar[i].fDegGalB<<" ";
						ATNFCheck<<pCfgPtr->pATNF->pPulsar[i].fP<<" ";
						ATNFCheck<<pCfgPtr->pATNF->pPulsar[i].fPdot<<" ";
						ATNFCheck<<pCfgPtr->pATNF->pPulsar[i].fRadioFlux1400<<" ";
						ATNFCheck<<pCfgPtr->pATNF->pPulsar[i].fDM<<" ";
	//					ATNFCheck<<pCfgPtr->pATNF->pPulsar[i].fDistDM<<" ";
						ATNFCheck<<pCfgPtr->pATNF->pPulsar[i].Name<<" ";

 						ATNFCheck<<endl;
					}

					pMaska[i]=true;
				}
			}
		}
	}

	if(pCfgPtr->Sim.bGaussianDensityMod)
	{
		//Clear the container
		for(IteratorMapTypeRadio Iterator = MapRadio.begin(); Iterator != MapRadio.end(); Iterator++)
		{
			Iterator->second=0.;
		}
		for(int i(0); i<pCfgPtr->pATNF->nTotalNum; i++)
		{
			if(pMaska[i])
			{
				for(IteratorMapTypeRadio Iterator = MapRadio.begin(); Iterator != MapRadio.end(); Iterator++)
				{
					Iterator->second += GaussianishDensityContribution(Iterator->first, ComputeSpacePointCoordiantesRadioSmooth( pCfgPtr->pATNF->pPulsar[i] ), pCfgPtr->Sim.fGaussianDensitySigma);//dodac tutaj sigmę
				}
			}
		}
	}



	//TODO 
	//map norm
	NormContainers();

	if(m_pCfgPtr->Sim.bComputeCatalogueMaxDist) fMaxDist = ComputeMaxDist(ID);
	else fMaxDist=m_pCfgPtr->Sim.fMaxDist;
	m_ComputeCM(nUsedATNFPulsars);

	if(pCfgPtr->Sim.bVerbose) 
	{
		cout<<"CMRadio: "<<CMRadio.fLogP<<" "<<CMRadio.fLogPdot<<" "<<CMRadio.fLogS1400<<endl;
		cout<<"CM4D: "<<CM4D.fLogP<<" "<<CM4D.fLogPdot<<" "<<CM4D.fLogS1400<<" "<<CM4D.fLogDM<<endl;
		cout<<"Used pulsars "<<nUsedATNFPulsars<<" Dead: "<<nDead<<endl;
		ATNFCheck.close();
		m_PrintUsedEntries(pMaska);
	}


	int nCells=0;
	for(IteratorMapTypeRadio Iterator = MapRadio.begin(); Iterator != MapRadio.end(); Iterator++) nCells++;

	if(pCfgPtr->Sim.bVerbose) 
	{
	cout<<"Effective catalogue size: "<<nUsedATNFPulsars<<endl;
	cout<<"Number od cells in RadioMap: "<<nCells<<endl;
    }

	delete [] pMaska;
}


void CObservations::NormContainers()
{
	Vector.Norm();

	//map norming
	RealType fTotal(0.);
	for(IteratorMapTypeRadio Iterator = MapRadio.begin(); Iterator != MapRadio.end(); Iterator++) 
	{
			fTotal+=Iterator->second;
	}
	for(IteratorMapTypeRadio Iterator = MapRadio.begin(); Iterator != MapRadio.end(); Iterator++) 
	{
			Iterator->second/=fTotal;
	}
	fTotal=0.;
	for(IteratorMapTypeSpace Iterator = MapSpace.begin(); Iterator != MapSpace.end(); Iterator++) 
	{
			fTotal+=Iterator->second;
	}
	for(IteratorMapTypeSpace Iterator = MapSpace.begin(); Iterator != MapSpace.end(); Iterator++) 
	{
			Iterator->second/=fTotal;
	}

	fTotal=0.;
	for(std::map<SProbabilitySpacePoint4D, RealType>::const_iterator Iterator = Map4D.begin(); Iterator!=Map4D.end(); Iterator++)
	{
			fTotal+=Iterator->second;
	}
	for(std::map<SProbabilitySpacePoint4D, RealType>::iterator Iterator = Map4D.begin(); Iterator!=Map4D.end(); Iterator++)
	{
			Iterator->second/=fTotal;
	}
}


void CObservations::m_PrintUsedEntries(const bool *pMaska)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	ofstream out("Used.ATNF.dat");
	//out<<"#P Pdot S1400"<<endl;
	SCatalogueEntry *pPSR = pCfgPtr->pATNF->pPulsar;
	for(int i(0);i<pCfgPtr->pATNF->nTotalNum; i++) if(pMaska[i])
	{
		out<<pPSR[i].fP<<" "<<pPSR[i].fPdot<<" "<<pPSR[i].fRadioFlux1400<<endl;
	}
	out.close();
}

void CObservations::m_AddToCM(const SCatalogueEntry Pulsar)
{
	CMRadio.fLogP += log10(Pulsar.fP);
	CMRadio.fLogPdot += log10(Pulsar.fPdot);
	CMRadio.fLogS1400 += log10(Pulsar.fRadioFlux1400);

	CM4D.fLogP += log10(Pulsar.fP);
	CM4D.fLogPdot += log10(Pulsar.fPdot);
	CM4D.fLogS1400 += log10(Pulsar.fRadioFlux1400);
	CM4D.fLogDM += log10(Pulsar.fDM);
}

void CObservations::m_ComputeCM(const int nTotalMass)
{
	CMRadio.fLogP /= nTotalMass;
	CMRadio.fLogPdot /= nTotalMass;
	CMRadio.fLogS1400 /= nTotalMass;

	CM4D.fLogP /= nTotalMass;
	CM4D.fLogPdot /= nTotalMass;
	CM4D.fLogS1400 /= nTotalMass;
	CM4D.fLogDM /= nTotalMass;

}

RealType CObservations::DistanceToCM(const CPulsar PSR)
{
	RealType DeltaLogP(0.);
	RealType DeltaLogPdot(0.);
	RealType DeltaLogS1400(0.);
	RealType DeltaLogDM(0.);
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density1x3DSparse  || m_pCfgPtr->Sim.ModelDensityContainer == Density1x3D || m_pCfgPtr->Sim.ModelDensityContainer == Density3x1D)
	{
		DeltaLogP = CMRadio.fLogP - log10(PSR.fPeriod);
		DeltaLogPdot = CMRadio.fLogPdot - log10(PSR.fPeriodDot);
		DeltaLogS1400 = CMRadio.fLogS1400 - log10(PSR.fS1400);
		//return sqrt(DeltaLogP*DeltaLogP + DeltaLogPdot*DeltaLogPdot + DeltaLogS1400*DeltaLogS1400) ;
		return (DeltaLogP*DeltaLogP + DeltaLogPdot*DeltaLogPdot + DeltaLogS1400*DeltaLogS1400) ;
		//return (DeltaLogP*DeltaLogP + DeltaLogPdot*DeltaLogPdot) ;
	}
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density4DSparse)
	{
		DeltaLogP = CM4D.fLogP - log10(PSR.fPeriod);
		DeltaLogPdot = CM4D.fLogPdot - log10(PSR.fPeriodDot);
		DeltaLogS1400 = CM4D.fLogS1400 - log10(PSR.fS1400);
		DeltaLogDM = CM4D.fLogDM - log10(PSR.fDM);
		//return sqrt(DeltaLogP*DeltaLogP + DeltaLogPdot*DeltaLogPdot + DeltaLogS1400*DeltaLogS1400 + DeltaLogDM*DeltaLogDM) ;
		return (DeltaLogP*DeltaLogP + DeltaLogPdot*DeltaLogPdot + DeltaLogS1400*DeltaLogS1400 + DeltaLogDM*DeltaLogDM) ;
		//return (DeltaLogP*DeltaLogP + DeltaLogPdot*DeltaLogPdot +  DeltaLogDM*DeltaLogDM) ;
	}

	//if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3DSparse || m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D )
	//TODO in this case an additional ful sum over distances to gal positions need to be compted
	return 0.;
}


void CObservations::m_AddToMap(const SProbabilitySpacePoint4D Point)
{
	pair< map<SProbabilitySpacePoint4D,RealType>::iterator,bool> RetVal;
	
	RetVal = Map4D.insert(pair<SProbabilitySpacePoint4D,int>(Point,1));
	if( !RetVal.second )
	{
		Map4D[Point]++;
		//if(Map4D[Point] > 20) ERROR("Current implementation do not support obs. greater then 20 in one volume element in the parameter space. Increse the resolution of the parameter space holder (radio).");
	}
}

void CObservations::m_AddToMap(const SProbabilitySpacePointRadio Point)
{
	pair< map<SProbabilitySpacePointRadio,RealType>::iterator,bool> RetVal;
	
	RetVal = MapRadio.insert(pair<SProbabilitySpacePointRadio,int>(Point,1));
	if( !RetVal.second )
	{
		MapRadio[Point]++;
		//if(MapRadio[Point] > 20) ERROR("Current implementation do not support obs. greater then 20 in one volume element in the parameter space. Increse the resolution of the parameter space holder (radio).");
	}
}

void CObservations::m_AddToMap(const SProbabilitySpacePointSpace Point)
{
	pair< map<SProbabilitySpacePointSpace,RealType>::iterator,bool> RetVal;
	
	RetVal = MapSpace.insert(pair<SProbabilitySpacePointSpace,int>(Point,1));
	if( !RetVal.second )
	{
		MapSpace[Point]++;
		//if(MapSpace[Point] > 20) ERROR("Current implementation do not support obs. greater then 20 in one volume element in the parameter space. Increse the resolution of the parameter space holder (spatial).");

	}
}


void  CObservations::Write(const char* Prefix)
{
	string sRadio="";
	if(Prefix!=NULL) sRadio += Prefix;
	sRadio += "Radio.Map.dat";	

	ofstream outRadio(sRadio.c_str());	
	RealType fLogP, fLogPdot, fLogS1400, fDegGalL, fKsi, fDM;
	for(IteratorMapTypeRadio Iterator = MapRadio.begin(); Iterator != MapRadio.end(); Iterator++) 
	{
		fLogP = m_pCfgPtr->Plot.fZeroLogP + Iterator->first.nLogP*(m_pCfgPtr->Plot.fMaxLogP-m_pCfgPtr->Plot.fZeroLogP)/m_pCfgPtr->Plot.nResLogP;
		fLogPdot = m_pCfgPtr->Plot.fZeroLogPdot + Iterator->first.nLogPdot*(m_pCfgPtr->Plot.fMaxLogPdot-m_pCfgPtr->Plot.fZeroLogPdot)/m_pCfgPtr->Plot.nResLogPdot;
		fLogS1400 = m_pCfgPtr->Plot.fZeroLogS1400 + Iterator->first.nLogS1400*(m_pCfgPtr->Plot.fMaxLogS1400-m_pCfgPtr->Plot.fZeroLogS1400)/m_pCfgPtr->Plot.nResLogS1400;

		outRadio<<fLogP<<" "<<fLogPdot<<" "<<fLogS1400<<" "<<Iterator->second<<endl;
	}
	outRadio.close();

	string sSpace = "";
	if(Prefix!=NULL) sSpace += Prefix;
	sSpace += "Space.Map.dat";
	
	ofstream outSpace(sSpace.c_str());
	for(IteratorMapTypeSpace Iterator = MapSpace.begin(); Iterator != MapSpace.end(); Iterator++) 
	{
		fDegGalL = Iterator->first.nGalL * 360. / m_pCfgPtr->Plot.nResGalL;
		fKsi = Iterator->first.nGalKsi * 2. / m_pCfgPtr->Plot.nResGalKsi - 1.;
		fDM = Iterator->first.nDM *  m_pCfgPtr->Plot.fMaxDM / m_pCfgPtr->Plot.nResDM;

		outSpace<<((fDegGalL>180.)?fDegGalL-360.:fDegGalL)<<" "<<fKsi<<" "<<fDM<<" "<<Iterator->second<<endl;

	}
	outSpace.close();
}

#include<Pulsar.hpp>
#include<RadioDetector.hpp>

void CTestObservations::Test()
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	int pksmb(0), req(0);

	CCumulativeDistribution CD_D(0.,450.,1000);
	CCumulativeDistribution CD_S(0.,80000.,10000);

	ofstream out("DM_S1400_diffproc.dat");
	CPulsar Psr;
	CRadioDetector Det;
	RealType dDMproc, dS1400proc;

	for(int i(0); i<pCfgPtr->pATNF->nTotalNum; i++)
	{
//		cout<<std::bitset<8>(pCfgPtr->pATNF->pPulsar[i].Surveys.Word & SurveyMask.Word)<<" "<<std::bitset<8>(SurveyMask.Word) <<endl;
//		cout<<std::bitset<16>(pCfgPtr->pATNF->pPulsar[i].DataFields.Word & pCfgPtr->Sim.DataFields.Word) << " "<<std::bitset<16>(pCfgPtr->Sim.DataFields.Word)<<endl;
		if((pCfgPtr->pATNF->pPulsar[i].Surveys.Word & SurveyMask.Word) == SurveyMask.Word )
		{
			pksmb++;
			if( (pCfgPtr->pATNF->pPulsar[i].DataFields.Word & pCfgPtr->Sim.DataFields.Word) == pCfgPtr->Sim.DataFields.Word)
			{
				req++;
				Psr.fPeriod = pCfgPtr->pATNF->pPulsar[i].fP;
				Psr.fPeriodDot = pCfgPtr->pATNF->pPulsar[i].fPdot;
				Psr.GalacticPos.fL = pCfgPtr->pATNF->pPulsar[i].fDegGalL*(M_PI/180.);
				Psr.GalacticPos.fB = pCfgPtr->pATNF->pPulsar[i].fDegGalB*(M_PI/180.);
				Psr.GalacticPos.fDist = pCfgPtr->pATNF->pPulsar[i].fDistDM;
				Psr.GalacticPosDeg.fL = pCfgPtr->pATNF->pPulsar[i].fDegGalL;
				Psr.GalacticPosDeg.fB = pCfgPtr->pATNF->pPulsar[i].fDegGalB;
				Psr.GalacticPosDeg.fDist = pCfgPtr->pATNF->pPulsar[i].fDistDM;
			
				Psr.fDM = Psr.DispersionMeasure(); 
				Psr.fL400= pow(10.,Psr.LogL400());//[mJy*kpc^2]

				Psr.fS1400 = Psr.RadioFlux(1400.);//mJy
				bool IsVisible=Det.CheckParkes(&Psr);
				dDMproc=100*((Psr.fDM - pCfgPtr->pATNF->pPulsar[i].fDM)/Psr.fDM);
				CD_D.AddPoint(dDMproc);
				dS1400proc= 100* (Psr.fS1400 - pCfgPtr->pATNF->pPulsar[i].fRadioFlux1400)/Psr.fS1400;
				CD_S.AddPoint(dS1400proc);
				//cout<<"DM: "<<Psr.fDM<<" "<<pCfgPtr->pATNF->pPulsar[i].fDM<<" "<<dDMproc<<endl;
				//cout<<"S1400: "<<Psr.fS1400<<" "<<pCfgPtr->pATNF->pPulsar[i].fRadioFlux1400*1e3<<" "<<dS1400proc<<endl;
				
				cout<<" "
				<<IsVisible<<" "
				<< pCfgPtr->pATNF->pPulsar[i].Name<<" "
				<<Psr.fPeriod<<"[s] "
				<<Psr.fPeriodDot<<"[s/s] "
				<<Psr.GalacticPos.fL<<"[rad] "
				<<Psr.GalacticPos.fB<<"[rad] "
				<<Psr.GalacticPos.fDist<<"[kpc] "
			//	<<Psr.fL400<<" "
				<<Psr.fS1400<<"[mJy] "
				<<fabs(dS1400proc)<<" [%] "
				<<Psr.fDM<<"[cm^-3 pc] "
				<<fabs(dDMproc)<<" [%] "
				<<endl;
				
				
				out<<fabs(dDMproc)<<" "<<fabs(dS1400proc)<<endl;

			}
		}
	}
	//cout<<pksmb<<" "<<req<<endl;
	CD_D.Cumulify();
	CD_D.Norm();
	CD_D.Write("Dyst.dDM.dat");

	CD_S.Cumulify();
	CD_S.Norm();
	CD_S.Write("Dyst.dS1400.dat");
	out.close();



}

RealType CObservations::ComputeMaxDist(const ESurveyID ID)
{
	nUsedATNFPulsars=0;
	SurveyID=ID;
	SurveyMask=SurveyID2Word(SurveyID);
	Survey.Init(SurveyID);
	CConfiguration *pCfgPtr=CConfiguration::GetInstance();

	m_pCfgPtr = CConfiguration::GetInstance();

	RealType fComputedDist, fComputedDM;
	RealType fMaxDist(0.), fMaxComputedDist(0.);
	//Not optimal but since all CObservations are constructed only once.
	//Count the number PSR's covered by given survey
	SCatalogueEntry PSR;
	if(m_pCfgPtr->Sim.bVerbose) cout<<"# L B DM ComputedDM Dist ComputedDist"<<endl;
	for(int i(0); i<pCfgPtr->pATNF->nTotalNum; i++)
	{
		if(m_pCfgPtr->Sim.bVerbose) 
		{
		//	cout<<"i: "<<i<<endl;
		//	cout<<std::bitset<8>(pCfgPtr->pATNF->pPulsar[i].Surveys.Word & SurveyMask.Word)<<" "<<std::bitset<8>(SurveyMask.Word) <<endl;
		}
		if((pCfgPtr->pATNF->pPulsar[i].Surveys.Word & SurveyMask.Word) == SurveyMask.Word )
		{
			if( (pCfgPtr->pATNF->pPulsar[i].DataFields.Word & pCfgPtr->Sim.DataFields.Word) == pCfgPtr->Sim.DataFields.Word)
			{
				nUsedATNFPulsars++;
				PSR = pCfgPtr->pATNF->pPulsar[i];

//				PSR.fDegGalB*(M_PI/180.);
//				PSR.fDegGalL*(M_PI/180.);
//				PSR.fDM;
//				PSR.fDist;
//				PSR.fDistDM;
				fComputedDist = fort_dmdsm(PSR.fDegGalL*(M_PI/180.),PSR.fDegGalB*(M_PI/180.),1,PSR.fDM,0.);
				fComputedDM = fort_dmdsm(PSR.fDegGalL*(M_PI/180.),PSR.fDegGalB*(M_PI/180.),-1,0.,PSR.fDist);
				if(m_pCfgPtr->Sim.bVerbose) 
				{
				cout<<PSR.fDegGalL<<" "
				    <<PSR.fDegGalB<<" "
                                    <<PSR.fDM<<" "
				    <<fComputedDM<<" "
                                    <<PSR.fDist<<" "
                                    <<fComputedDist<<endl;
				}
				if(PSR.fDist>fMaxDist) fMaxDist = PSR.fDist;
				if(fComputedDist>fMaxComputedDist) fMaxComputedDist=fComputedDist;
			}
		}
	}
	if(m_pCfgPtr->Sim.bVerbose) 
	{
		cout<<"# All:"<<pCfgPtr->pATNF->nTotalNum<<endl;
		cout<<"# Used: "<<nUsedATNFPulsars<<endl;
		cout<<"# MaxDist: "<<fMaxDist<<endl;
		cout<<"# MaxComputedDist: "<<fMaxComputedDist<<endl;
		cout<<"# Result: "<<(fMaxDist > fMaxComputedDist ? fMaxDist : fMaxComputedDist)<<endl;
	}
	return fMaxDist > fMaxComputedDist ? fMaxDist : fMaxComputedDist ;
	//cout<<"Used pulsars "<<nUsedATNFPulsars<<endl;
	
}


