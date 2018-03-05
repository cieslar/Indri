#include<MCMC.hpp>
#include<MyMath.hpp>

#include<Randomizer.hpp>
#include<Statistics.hpp>

#include<fstream>
#include<iostream>

#include<cmath>
#include<string>
#include<sstream>

#include <H5Cpp.h>

using namespace std;

void CModel::PrintPulsar(const int nNum) const
{
	m_pPop->pPulsar[nNum].PrintDescriptiveStarInfo(std::cout);
}


void CModel::Test()
{
//	m_pCfgPtr->pATNF->Print();

	m_ppObs[0]->Write();

}

CModel::CModel()
{

	m_pCfgPtr = CConfiguration::GetInstance();
	m_pPop = new CPopulation();
	if(m_pCfgPtr->Sim.bVerbose) cout<<"Created population holder"<<endl;
	m_nObsNum = m_pCfgPtr->SurveySelection.nNum;
	m_pCfgPtr->SurveysCuts.pMaxSurveyDist = new RealType[m_nObsNum];
	m_pCfgPtr->SurveysCuts.bInitiated = true;

	m_pRadioDet = new CRadioDetector();
	if(m_pCfgPtr->Sim.bVerbose) cout<<"Detector constructed"<<endl;

	m_ppObs = new CObservations*[m_nObsNum];
	m_ppDensity= new CDensity*[m_nObsNum];

	//cout<< m_pCfgPtr->SurveySelection.nNum<<endl;

	for(int i(0); i<m_nObsNum; i++)
	{
		m_ppObs[i] = new CObservations(m_pCfgPtr->SurveySelection.pID[i]);
		if(m_pCfgPtr->Sim.bVerbose) cout<<"Created "<<i<<"-observation holder"<<endl;
		m_pCfgPtr->SurveysCuts.pMaxSurveyDist[i] = m_ppObs[i]->fMaxDist;
		if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D || m_pCfgPtr->Sim.ModelDensityContainer == Density1x3D )
		{
			m_ppDensity[i] = dynamic_cast<CDensity*>(new CDensity2x3D());
			if(m_pCfgPtr->Sim.bVerbose) cout<<"Created models' "<<i<<"-density holder"<<endl;
		}

		if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3DSparse || m_pCfgPtr->Sim.ModelDensityContainer == Density1x3DSparse )
		{
			m_ppDensity[i] = dynamic_cast<CDensity*>(new CDensity2x3DSparse(m_ppObs[i] ));
			if(m_pCfgPtr->Sim.bVerbose) cout<<"Created models' "<<i<<"-density holder"<<endl;
		}

		if(m_pCfgPtr->Sim.ModelDensityContainer == Density4DSparse )
		{
			m_ppDensity[i] = dynamic_cast<CDensity*>(new CDensity4DSparse(m_ppObs[i] ));
			if(m_pCfgPtr->Sim.bVerbose) cout<<"Created models' "<<i<<"-density holder"<<endl;
		}
		
		if(m_pCfgPtr->Sim.ModelDensityContainer == Density3x1D )
		{
			m_ppDensity[i] = dynamic_cast<CDensity*>(new CDensity3x1D());
			if(m_pCfgPtr->Sim.bVerbose) cout<<"Created models' "<<i<<"-density holder"<<endl;
		}

	}

	//Zero();//is it needed?
	if(m_pCfgPtr->Sim.ExtremumMethod == KolmogorovSmirnov) m_fWorstCaseLnLikelihood = 1.;
	m_fWorstCaseLnLikelihood = LnLikelihood();
	fTotalPenalty=0.;
	//ERROR("DUPA");
}

void CModel::Check()
{
	for(int i(0); i<m_nObsNum; i++)
	{
		m_ppDensity[i]->Check();
	}


}
CModel::~CModel()
{
	delete m_pRadioDet;
	delete m_pPop;
	for(int i(0);i<m_nObsNum;i++)
	{
		delete m_ppObs[i];
		delete m_ppDensity[i];
	}
	delete [] m_ppObs;
	delete [] m_ppDensity;

}


void CModel::FullCycle(const char* HDFFileName)
{
	for(int i(0); i<m_nObsNum; i++)
	{
		m_ppDensity[i]->Zero();
	}


	for(int i(0);i<m_pCfgPtr->Sim.nNumberOfIteration;i++)
	{
		if(m_pCfgPtr->Sim.bVerbose)
		{
			std::cout<<"Model's iteration: "<<i<<std::endl;
		}
		m_pPop->Zero();
		m_pPop->Evolve();
		FillDensityContainers();
		if(HDFFileName != NULL) m_pPop->WriteHDF5(HDFFileName,i);

		m_pCfgPtr->Sim.bRecycleDynamics=true;
	}

	NormalizeDensityContainers();
	ComputeMinimizationValue();
}



void CModel::ZeroPop()
{
	m_pPop->Zero();
}

void CModel::Zero()
{
	m_pPop->Zero();
	for(int i(0); i<m_nObsNum; i++)
	{
		m_ppDensity[i]->Zero();
	}

}

void CModel::Evolve()
{
	m_pPop->Evolve();//FIXME upewnić się czy to wszystko oblicza
}



void CModel::FillDensityContainers()
{
    int nNotVisible(0);
    int nLowB(0);
    int nLost(0);
    int nCount(0);
    int nN(0);

	for(int k(0);k<m_nObsNum;k++)
	{
        nNotVisible=nLowB=nLost=nCount=nN=0;

		for(int i(0);i<m_pCfgPtr->Sim.nMaxPulsarNum; i++)
		{
			if( !m_pPop->pPulsar[i].bIsLost) 
			{
				//same condition as the preselction process from ATNF Catalogue - to avoid the MSP pulsars
				if( m_pPop->pPulsar[i].fBField >= 1e10 )
                {
				//whis may be a redundant check
				//TODO clean in the future
				//TODO TODO this generates additional traffic in the Survey.init()?! FIX IT!
                    if( m_pRadioDet->IsVisible(&(m_pPop->pPulsar[i]),m_ppObs[k]->SurveyID) )
                    {
                        nCount++;
                        if(m_pCfgPtr->Sim.bGaussianDensityMod)
                        {
                            for(IteratorMapTypeRadio Iterator = m_ppObs[k]->MapRadio.begin(); Iterator != m_ppObs[k]->MapRadio.end(); Iterator++) 
                            {
                                m_ppDensity[k]->AddValueAtPoint( Iterator->first, 
                                                                 GaussianishDensityContribution(Iterator->first, 
                                                                                                ComputeSpacePointCoordiantesRadioSmooth( &(m_pPop->pPulsar[i]) ), 
                                                                                                m_pCfgPtr->Sim.fGaussianDensitySigma) );

                                nN++;
                            }
                        }
                        else
                        {
                            m_ppDensity[k]->AddPoint(m_pPop->pPulsar[i]);
                        }
                    
                    }
                    else
                    {
                        nNotVisible++;
                    }
                }
                else
                {
                    nLowB++;
                }
			}
            else
            {
                nLost++;
            }
		}
        if(m_pCfgPtr->Sim.bVerbose)
        {
            cout<<"Num. of pulsars used for the lnL computations: "<<nCount<<" "<<nN/nCount<<endl;
            cout<<"Not visible: "<<100*static_cast<float>(nNotVisible)/m_pCfgPtr->Sim.nMaxPulsarNum<<"%"<<endl;
            cout<<"Too low B: "<<100*static_cast<float>(nLowB)/m_pCfgPtr->Sim.nMaxPulsarNum<<"%"<<endl;
            cout<<"Lost: "<<100*static_cast<float>(nLost)/m_pCfgPtr->Sim.nMaxPulsarNum <<"%"<<endl;
        }
	}

}

void CModel::WriteAsCatalogue(const char* FileName, const int nNumber)
{
	WARNING("The artificial catalogue is set to the Parkes Multibeam Survey.");
	//Need to correct the AsCatalogueEntry() to get as input the SurveyID
	//Set to Parkes for now
	
	ofstream OutFile(FileName);
	int nWritten(0);
	for(int k(0);k<m_nObsNum;k++)//will be incorrect if more then Parkes
	{
		for(int i(0);i<m_pCfgPtr->Sim.nMaxPulsarNum && nWritten < nNumber; i++)
		{
			if( !m_pPop->pPulsar[i].bIsLost) 
			{
				//same condition as the preselction process from ATNF Catalogue - to avoid the MSP pulsars
				if( m_pPop->pPulsar[i].fBField >= 1e10 )
				//whis may be a redundant check
				//TODO clean in the future
				if( m_pRadioDet->IsVisible(&(m_pPop->pPulsar[i]),m_ppObs[k]->SurveyID) )
				{
					OutFile<<m_pPop->pPulsar[i].AsCatalogueEntry(i);
					nWritten++;
				}
			}
		}
	}
	OutFile.close();
}


void CModel::NormalizeDensityContainers()
{
	for(int k(0);k<m_nObsNum;k++)
	{
		if(m_pCfgPtr->Sim.bPhysicalDensity) m_ppDensity[k]->Densify();
		m_ppDensity[k]->Norm();
        //for(IteratorMapTypeRadio Iterator = m_ppObs[k]->MapRadio.begin(); Iterator != m_ppObs[k]->MapRadio.end(); Iterator++) {
        //    cout<<"\t"<<m_ppDensity[k]->GetDensityRadio(Iterator->first)<<" "<<Iterator->second;
        //}
	}
}

//reworked
void CModel::Densify()
{
	RealType Eff(0.), fPenaltyOk, fPenaltyRest;
	int nOK(0);
	fTotalPenalty=0.;
	for(int k(0);k<m_nObsNum;k++)
	{
		nOK=0;
		fPenaltyOk=fPenaltyRest=0.;
		for(int i(0);i<m_pCfgPtr->Sim.nMaxPulsarNum; i++)
		{

			if( !m_pPop->pPulsar[i].bIsLost) 
			{
				//TODO czy to nie jest tozsame z bIsLost??
				if( m_pRadioDet->IsVisible(&(m_pPop->pPulsar[i]),m_ppObs[k]->SurveyID) )
				{
					nOK++;
					m_ppDensity[k]->AddPoint(m_pPop->pPulsar[i]);
					fPenaltyOk += m_ppObs[k]->DistanceToCM(m_pPop->pPulsar[i]);
				}
				else
				{
					fPenaltyRest += m_ppObs[k]->DistanceToCM(m_pPop->pPulsar[i]);
				}
			}
		}
		if(m_pCfgPtr->Sim.bPhysicalDensity) m_ppDensity[k]->Densify();
		m_ppDensity[k]->Norm();
		if(nOK>0)
		{
			Eff+=m_ppDensity[k]->GetEfficiency()/m_nObsNum;
			fPenaltyOk /= nOK;
		}
		fPenaltyRest/=(m_pCfgPtr->Sim.nMaxPulsarNum-nOK);
	
		fTotalPenalty += m_pCfgPtr->Sim.fLnLikelihoodPenaltyCoefficient * fPenaltyOk; 
		if(nOK==0) fTotalPenalty += (m_pCfgPtr->Sim.fLnLikelihoodPenaltyCoefficient) * fPenaltyRest;
	}
	Efficiency=Eff;	
}

void CModel::Write(const char* Prefix)
{	
	string Pre;
	if(Prefix==NULL) Pre="";
	else Pre=Prefix;

	for(int k(0);k<m_nObsNum;k++)
	{
		
		string Name("Survey");
		Name = Prefix + Name + static_cast<char>(k+48);
		m_ppDensity[k]->Write(Name.c_str());
	}
	string Pop=Prefix;
	Pop+="pop.dat";
	m_pPop->Write(Prefix);

	string sCat;
	for(int i(0);i<m_nObsNum;i++)
	{
		sCat=Prefix;
		sCat+="-";
		sCat+=std::to_string(i);
		m_ppObs[i]->Write(sCat.c_str());
	}
}

void CModel::ComputeMinimizationValue()
{
	m_pCfgPtr->Model.fLnLikelihood = LnLikelihood();
	m_pCfgPtr->Model.fDValue = DValue();
}

RealType CModel::DValue()
{
	RealType fDValue(0.);
	if( m_pCfgPtr->Sim.ModelDensityContainer == Density3x1D )
	{
		for(int i(0);i<m_nObsNum;i++)
		{
			fDValue += DValue3D(&(static_cast<CDensity3x1D*>(m_ppDensity[i])->Vector),&(m_ppObs[i]->Vector));
		}
	}
	else
	{
		//WARNING("Density container not compatible. DValue defaulted to 1.");
		fDValue=1.;
	}
	//return log10(fDValue);
	return fDValue;
}


//RealType CModel::LnLikelihoodFunction

RealType CModel::LnLikelihood()
{
	RealType fLnLikelihood(0.);
	RealType fDensity(0.);
	int nNum;

	bool bWorstCase(false);

	RealType fGaussSigma(m_pCfgPtr->Sim.fLnLikelihoodSigma);

	//cout<<"LnLikelihood Model.DensityContainer: "<<m_pCfgPtr->Sim.ModelDensityContainer<<endl;

	if(m_pCfgPtr->Sim.ModelDensityContainer == Density3x1D )
	for(int i(0);i<m_nObsNum;i++)
	{
		for(int j(0);j< m_pCfgPtr->Plot.nResLogP;j++)
		{
			fLnLikelihood += LnGaussProbability( static_cast<CDensity3x1D*>(m_ppDensity[i])->Vector.pPeriodLog[j], m_ppObs[i]->Vector.pPeriodLog[j],fGaussSigma);
		}

		for(int j(0);j< m_pCfgPtr->Plot.nResLogPdot;j++)
		{
			fLnLikelihood += LnGaussProbability( static_cast<CDensity3x1D*>(m_ppDensity[i])->Vector.pPeriodDotLog[j], m_ppObs[i]->Vector.pPeriodDotLog[j],fGaussSigma );

		}
		
		for(int j(0);j< m_pCfgPtr->Plot.nResLogS1400;j++)
		{
			fLnLikelihood += LnGaussProbability( static_cast<CDensity3x1D*>(m_ppDensity[i])->Vector.pS1400Log[j], m_ppObs[i]->Vector.pS1400Log[j],fGaussSigma );
		}


	}

	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D || m_pCfgPtr->Sim.ModelDensityContainer == Density1x3D )
	for(int i(0);i<m_nObsNum;i++)
	{
		for(IteratorMapTypeRadio Iterator = m_ppObs[i]->MapRadio.begin(); Iterator != m_ppObs[i]->MapRadio.end(); Iterator++) 
		{
			if(Iterator->second > 0)
			{
				fLnLikelihood += LnGaussProbability(  m_ppDensity[i]->GetDensityRadio(Iterator->first) ,Iterator->second,fGaussSigma );
			}
            else
            {
                if(m_pCfgPtr->Sim.bVerbose)
			    cout<<"("<<Iterator->first.nLogP<<","<<Iterator->first.nLogPdot<<","<<Iterator->first.nLogS1400<<") "<<Iterator->second<<" "<<m_ppDensity[i]->GetDensityRadio(Iterator->first) <<endl;
            }
		}
		if( m_pCfgPtr->Sim.ModelDensityContainer == Density2x3D )
		for(IteratorMapTypeSpace Iterator = m_ppObs[i]->MapSpace.begin(); Iterator != m_ppObs[i]->MapSpace.end(); Iterator++) 
		{
			if(Iterator->second > 0)
			{
				fLnLikelihood += LnGaussProbability(  m_ppDensity[i]->GetDensitySpace(Iterator->first) ,Iterator->second,fGaussSigma );
			}
		}
	}
	//for each survey compute the Likelihood and return
	//I will use log of the Likelihood
	//cout<<fLnLikelihood<<endl;
		//check if the Density is zero - use the fill value;
	if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3DSparse || m_pCfgPtr->Sim.ModelDensityContainer == Density1x3DSparse )
	for(int i(0);i<m_nObsNum;i++)
	{
		std::map<SProbabilitySpacePointRadio, RealType>::const_iterator ObservationsIterator;
		std::map<SProbabilitySpacePointRadio, RealType>::const_iterator ModelIterator;

		ObservationsIterator = m_ppObs[i]->MapRadio.begin();
		CDensity2x3DSparse* pDensity = dynamic_cast<CDensity2x3DSparse*>(m_ppDensity[i]);
		ModelIterator        = pDensity->MapRadio.begin();

		//They have the same keys in the same order
		while(ObservationsIterator!=m_ppObs[i]->MapRadio.end())
		{
			fLnLikelihood += LnGaussProbability( ModelIterator->second, ObservationsIterator->second,fGaussSigma );

			++ObservationsIterator;
			++ModelIterator;
		}
		if(m_pCfgPtr->Sim.ModelDensityContainer == Density2x3DSparse)
		{
			std::map<SProbabilitySpacePointSpace, RealType>::const_iterator ObservationsIterator;
			std::map<SProbabilitySpacePointSpace, RealType>::const_iterator ModelIterator;

			ObservationsIterator = m_ppObs[i]->MapSpace.begin();
			CDensity2x3DSparse* pDensity = dynamic_cast<CDensity2x3DSparse*>(m_ppDensity[i]);
			ModelIterator        = pDensity->MapSpace.begin();

			//They have the same keys in the same order
			while(ObservationsIterator!=m_ppObs[i]->MapSpace.end())
			{
				fLnLikelihood += LnGaussProbability( ModelIterator->second, ObservationsIterator->second,fGaussSigma );
			
				++ObservationsIterator;
				++ModelIterator;
			}
		}
	}

	if(m_pCfgPtr->Sim.ModelDensityContainer == Density4DSparse)
	for(int i(0);i<m_nObsNum;i++)
	{
		std::map<SProbabilitySpacePoint4D, RealType>::const_iterator ObservationsIterator;
		std::map<SProbabilitySpacePoint4D, RealType>::const_iterator ModelIterator;

		ObservationsIterator = m_ppObs[i]->Map4D.begin();
		CDensity4DSparse* pDensity = dynamic_cast<CDensity4DSparse*>(m_ppDensity[i]);
		ModelIterator        = pDensity->Map4D.begin();

		//They have the same keys in the same order
		while(ObservationsIterator!=m_ppObs[i]->Map4D.end())
		{
			fLnLikelihood += LnGaussProbability( ModelIterator->second, ObservationsIterator->second,fGaussSigma );
			++ObservationsIterator;
			++ModelIterator;
		}

	}
	if(m_pCfgPtr->Sim.bLnLikelihoodPenalty) 
	{
		fLnLikelihood *= fTotalPenalty;
	}

	if(m_pCfgPtr->Sim.bStrictLikelihood && bWorstCase) fLnLikelihood=m_fWorstCaseLnLikelihood;
	return fLnLikelihood;
}

void CModel::ReadBinary(const char* PopFile)
{
	m_pPop->ReadBinary(PopFile);

}

void CModel::WriteBinaryHDF5(const char* FileName, const int nIteration)
{
	m_pPop->WriteHDF5(FileName,nIteration);
}
void CModel::WritePartialBinaryHDF5(const char* FileName)
{
	//m_pPop->WriteHDF5(FileName);
}
//void CModel::ReadBinaryHDF5(const char* FileName)
//{
//	m_pPop->ReadBinaryHDF5(FileName);
//}




void CTestModel::Test()
{
	Zero();
	Evolve();
	Densify();
//	Write();
//	CModel::Test();
	//cout<<m_pCfgPtr->Model.fBDecayTimeDeltaMean<<" "<<LnLikelihood()<<endl;

}








void CMCMC::PopReadBinary(const char* PopFile)
{
	m_pModel->ReadBinary(PopFile);
}


CMCMC::CMCMC()
{
	m_pCfgPtr = CConfiguration::GetInstance();
	m_pRanPtr = m_pCfgPtr->pRan;

	m_nNumOfChains=m_pCfgPtr->Sim.nNumOfChains;
	m_nNumPerChain=m_pCfgPtr->Sim.nNumPerChain;
	m_bChainMemOn=false;
	m_pModel = new CModel();
	
	m_CreateEmptyChains();

}

void CMCMC::m_CreateEmptyChains()
{
	//utworzenie lancuchow i wyzerowanie
	m_pChain = new SModelParameters*[m_pCfgPtr->Sim.nNumOfChains];
	for(int i(0); i<m_pCfgPtr->Sim.nNumOfChains; i++)
	{
		m_pChain[i] = new SModelParameters[m_pCfgPtr->Sim.nNumPerChain];
		for(int j(0); j<m_pCfgPtr->Sim.nNumPerChain; j++)
		{
			//Zero(m_pChain[i][j]);//Dopisać konstruktor zerujący
		}
		
	}
	m_bChainMemOn=true;

}


void CMCMC::m_DeleteChains()
{
	if(m_bChainMemOn)
	{
		for(int i(0); i<m_pCfgPtr->Sim.nNumOfChains; i++)
		{
			delete [] m_pChain[i];
		}
		delete [] m_pChain;	
		m_bChainMemOn=false;
	}
}


CMCMC::~CMCMC()
{
	m_DeleteChains();	
}







SModelParameters CMCMC::m_StartChain()
{
	///Draw from the priror distribution
    SModelParameters CurrentParams;
    if(m_pCfgPtr->Sim.bStartMCMCFromModel)
    {
        if(m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Starting MCMC from the Model"<<endl;
        CurrentParams = m_pCfgPtr->Model; 
    }
    else
    {
        CurrentParams = InitialDistribution();
    }

	m_ComputeValue(CurrentParams);
	
	//cout<<CurrentParams.fLnLikelihood<<" "<<CurrentParams.fBDecayTimeDeltaMean<<" "<< CurrentParams.fBDecayTimeDeltaSigma<<endl;
	return CurrentParams;
}

SModelParameters CMCMC::m_MakeLink(SModelParameters Previous)
{
	SModelParameters Next = SamplingDistribution(Previous);

	//Added due to the intrinsic likelihood uncertainty 
	//m_ComputeValue(Previous);

	m_ComputeValue(Next);

	if(m_pCfgPtr->Sim.ExtremumMethod == Likelihood )
	{
		RealType fLnR = Next.fLnLikelihood -  Previous.fLnLikelihood;
		RealType fLnK = log(m_pRanPtr->FlatLeftOpen());

        if(m_pCfgPtr->Sim.bMCMCVerbose) cout<<(fLnR >= 0.)<<" "<<(fLnK <= fLnR)<<" "<<Previous.fLnLikelihood<<" "<<Next.fLnLikelihood<<endl;

		if( (fLnR >= 0.) || (fLnK <= fLnR) )
		{
			return Next;
		}
		return Previous;
	}

	if(m_pCfgPtr->Sim.ExtremumMethod == KolmogorovSmirnov )
	{
		if(Next.fDValue < Previous.fDValue)
		{
			return Next;
		}
		return Previous;
	}
	
	WARNING("No valid extremum method!");
	return Previous;
}


void CMCMC::MakeChains()
{
	int nStepResizingPoint(m_pCfgPtr->Sim.nScaleEveryJumps);

	for(int i(0); i<m_nNumOfChains; i++)
	{
		if(m_pCfgPtr->Sim.bMCMCVerbose)
		{
			cout<<"Chain: "<<i<<" Link: "<<0<<endl;
			cout.flush();
		}
		m_pChain[i][0] = m_StartChain();
		for(int j(1); j<m_pCfgPtr->Sim.nNumPerChain; j++)
		{
			if(m_pCfgPtr->Sim.bMCMCVerbose)
			{
				cout<<"Chain: "<<i<<" Link: "<<j<<endl;
				cout.flush();
			}

            if(m_pCfgPtr->Sim.bScaleJumps)
			if( (j/nStepResizingPoint)*nStepResizingPoint == j)
			{
				ScaleJumpSigmas();

				if(m_pCfgPtr->Sim.bMCMCVerbose)
				{
					cout<<"Jump simgas scalled"<<endl;
					cout.flush();
				}
			}
	
			m_pChain[i][j] = m_MakeLink( m_pChain[i][j-1] );
		}
	}
}



void CMCMC::Test()
{
	//MakeChains();
	//m_StartChain();

	SModelParameters Params;
	ofstream out("InitParamsTest.txt");
	for(int i(0);i<1000000;++i)
	{
		Params=InitialDistribution();
		out<<Params.fLPowerLawGamma<<" ";
		out<<Params.fLPowerLawAlpha<<" ";
		out<<Params.fLPowerLawBeta<<" ";
		out<<Params.fBDecayDelta<<" ";
		out<<Params.fBDecayKappa<<" ";
		out<<endl;
	}
	out.close();

	out.open("PropParamsTest.txt");
	for(int i(0);i<1000000;++i)
	{
		Params=SamplingDistribution(Params);
		out<<Params.fLPowerLawGamma<<" ";
		out<<Params.fLPowerLawAlpha<<" ";
		out<<Params.fLPowerLawBeta<<" ";
		out<<Params.fBDecayDelta<<" ";
		out<<Params.fBDecayKappa<<" ";
		out<<endl;


	}
	out.close();
}
//#include<Statistics.hpp>
void CMCMC::MakeDensitySpace()
{
	//CAreaDensity AD(m_pCfgPtr->ModelConstraints.Min.fBDecayTimeDeltaMean, m_pCfgPtr->ModelConstraints.Min.fBDecayTimeDeltaSigma,  m_pCfgPtr->ModelConstraints.Max.fBDecayTimeDeltaMean, m_pCfgPtr->ModelConstraints.Max.fBDecayTimeDeltaSigma, 128, 128);



	//for(int i(0); i<m_nNumOfChains; i++) AD.AddPoint(m_pChain[i][m_pCfgPtr->Sim.nNumPerChain-1].fBDecayTimeDeltaMean, m_pChain[i][m_pCfgPtr->Sim.nNumPerChain-1].fBDecayTimeDeltaSigma);

//	AD.Densify();
//	AD.Norm();
//	AD.Write("ParameterSpace.mtx");

}

void CMCMC::WriteWholeChains()
{
		



}

void CMCMC::Write()
{
	Write(m_pCfgPtr->Sim.nNumPerChain);
}

void CMCMC::Write(const int nPointOfChain)
{
	ofstream out("Results/LastChain.dat");
	int nNum = MINIMUM(nPointOfChain-1,m_pCfgPtr->Sim.nNumPerChain-1);
	
//	CHistogram *HistBDecayTimeDeltaMean;
//	CHistogram *HistBDecayTimeDeltaSigma;

	

	SModelParameters Avg;
	SModelParameters Sgm;

	SModelParameters LogAvg;
	SModelParameters LogSgm;

	RealType CovGammaAlpha=0.;
	RealType CovLogGammaAlpha=0.;
	RealType CovGammaBeta=0.;
	RealType CovLogGammaBeta=0.;
	RealType CovGammaKappa=0.;
	RealType CovLogGammaKappa=0.;
	RealType CovGammaDelta=0.;
	RealType CovGammaLogDelta=0.;
	RealType CovLogGammaDelta=0.;
	RealType CovLogGammaLogDelta=0.;
	RealType CovAlphaBeta=0.;
	RealType CovAlphaKappa=0.;
	RealType CovAlphaDelta=0.;
	RealType CovAlphaLogDelta=0.;
	RealType CovBetaKappa=0.;
	RealType CovBetaDelta=0.;
	RealType CovBetaLogDelta=0.;
	RealType CovKappaDelta=0.;
	RealType CovKappaLogDelta=0.;
	//if (m_pCfgPtr->ModelConstraints.bBDecayTimeDeltaMean) HistBDecayTimeDeltaMean = new CHistogram(m_pCfgPtr->ModelConstraints.Min.fBDecayTimeDeltaMean, m_pCfgPtr->ModelConstraints.Max.fBDecayTimeDeltaMean, 20);
	//if (m_pCfgPtr->ModelConstraints.bBDecayTimeDeltaSigma) HistBDecayTimeDeltaSigma = new CHistogram(m_pCfgPtr->ModelConstraints.Min.fBDecayTimeDeltaSigma, m_pCfgPtr->ModelConstraints.Max.fBDecayTimeDeltaSigma, 256);

	//The header
//	out<<"# Chains: "<<m_pCfgPtr->Sim.nNumOfChains<<" Links: "<<m_pCfgPtr->Sim.nNumPerChain<<" OutputLink: "<<nNum<<endl;
//	out<<"# ";
//	if (m_pCfgPtr->ModelConstraints.bBDecayTimeDeltaMean)  out<<"BDecayTimeDeltaMean[Myr] ";
//	if (m_pCfgPtr->ModelConstraints.bBDecayTimeDeltaSigma) out<<"BDecayTimeDeltaSigma[Myr] ";
//	out<<"LnLikelihood"<<endl;
//
	for(int i(0); i<m_pCfgPtr->Sim.nNumOfChains; i++)
	{
		Avg = Avg + m_pChain[i][nNum];

		LogAvg = LogAvg + log10(m_pChain[i][nNum]);

		out<<m_pChain[i][nNum].fLPowerLawGamma<<" ";
		out<<m_pChain[i][nNum].fLPowerLawAlpha<<" ";
		out<<m_pChain[i][nNum].fLPowerLawBeta<<" ";
		out<<m_pChain[i][nNum].fBDecayDelta<<" ";
		out<<m_pChain[i][nNum].fBDecayKappa<<" ";
		out<<m_pChain[i][nNum].fLnLikelihood<<" ";
		out<<endl;
	}	
	out.close();
	
	Avg = Avg * (1./m_pCfgPtr->Sim.nNumOfChains);
	LogAvg = LogAvg * (1./m_pCfgPtr->Sim.nNumOfChains);
	for(int i(0); i<m_pCfgPtr->Sim.nNumOfChains; i++)
	{
		Sgm = Sgm + pow((m_pChain[i][nNum] - Avg),2.);
		LogSgm = LogSgm + pow((log10(m_pChain[i][nNum]) - LogAvg),2.);

		CovGammaAlpha       += (m_pChain[i][nNum].fLPowerLawGamma - Avg.fLPowerLawGamma)*(m_pChain[i][nNum].fLPowerLawAlpha - Avg.fLPowerLawAlpha);
		CovLogGammaAlpha    += (log10(m_pChain[i][nNum].fLPowerLawGamma) - LogAvg.fLPowerLawGamma)*(m_pChain[i][nNum].fLPowerLawAlpha - Avg.fLPowerLawAlpha);
		CovGammaBeta        += (m_pChain[i][nNum].fLPowerLawGamma - Avg.fLPowerLawGamma)*(m_pChain[i][nNum].fLPowerLawBeta - Avg.fLPowerLawBeta);
		CovLogGammaBeta     += (log10(m_pChain[i][nNum].fLPowerLawGamma) - LogAvg.fLPowerLawGamma)*(m_pChain[i][nNum].fLPowerLawBeta - Avg.fLPowerLawBeta);
		CovGammaKappa       += (m_pChain[i][nNum].fLPowerLawGamma - Avg.fLPowerLawGamma)*(m_pChain[i][nNum].fBDecayKappa - Avg.fBDecayKappa);
		CovLogGammaKappa    += (log10(m_pChain[i][nNum].fLPowerLawGamma) - LogAvg.fLPowerLawGamma)*(m_pChain[i][nNum].fBDecayKappa - Avg.fBDecayKappa);
		CovGammaDelta       += (m_pChain[i][nNum].fLPowerLawGamma - Avg.fLPowerLawGamma)*(m_pChain[i][nNum].fBDecayDelta - Avg.fBDecayDelta);
		CovGammaLogDelta    += (m_pChain[i][nNum].fLPowerLawGamma - Avg.fLPowerLawGamma)*(log10(m_pChain[i][nNum].fBDecayDelta) - LogAvg.fBDecayDelta);
		CovLogGammaDelta    += (log10(m_pChain[i][nNum].fLPowerLawGamma) - LogAvg.fLPowerLawGamma)*(m_pChain[i][nNum].fBDecayDelta - Avg.fBDecayDelta);
		CovLogGammaLogDelta += (log10(m_pChain[i][nNum].fLPowerLawGamma) - LogAvg.fLPowerLawGamma)*(log10(m_pChain[i][nNum].fBDecayDelta) - LogAvg.fBDecayDelta);
		
		CovAlphaBeta        += (m_pChain[i][nNum].fLPowerLawAlpha - Avg.fLPowerLawAlpha)*(m_pChain[i][nNum].fLPowerLawBeta - Avg.fLPowerLawBeta);
		CovAlphaKappa       += (m_pChain[i][nNum].fLPowerLawAlpha - Avg.fLPowerLawAlpha)*(m_pChain[i][nNum].fBDecayKappa - Avg.fBDecayKappa);
		CovAlphaDelta       += (m_pChain[i][nNum].fLPowerLawAlpha - Avg.fLPowerLawAlpha)*(m_pChain[i][nNum].fBDecayDelta - Avg.fBDecayDelta);
		CovAlphaLogDelta    += (m_pChain[i][nNum].fLPowerLawAlpha - Avg.fLPowerLawAlpha)*(log10(m_pChain[i][nNum].fBDecayDelta) - LogAvg.fBDecayDelta);

		CovBetaKappa        += (m_pChain[i][nNum].fLPowerLawBeta - Avg.fLPowerLawBeta)*(m_pChain[i][nNum].fBDecayKappa - Avg.fBDecayKappa);
		CovBetaDelta        += (m_pChain[i][nNum].fLPowerLawBeta - Avg.fLPowerLawBeta)*(m_pChain[i][nNum].fBDecayDelta - Avg.fBDecayDelta);
		CovBetaLogDelta     += (m_pChain[i][nNum].fLPowerLawBeta - Avg.fLPowerLawBeta)*(log10(m_pChain[i][nNum].fBDecayDelta) - LogAvg.fBDecayDelta);

		CovKappaDelta        += (m_pChain[i][nNum].fBDecayKappa - Avg.fBDecayKappa)*(m_pChain[i][nNum].fBDecayDelta - Avg.fBDecayDelta);
		CovKappaLogDelta     += (m_pChain[i][nNum].fBDecayKappa - Avg.fBDecayKappa)*(log10(m_pChain[i][nNum].fBDecayDelta) - LogAvg.fBDecayDelta);

	}

	RealType fN(m_pCfgPtr->Sim.nNumOfChains);
	
	Sgm = sqrt(Sgm * (1./fN) );
	LogSgm = sqrt(LogSgm * (1./fN) );


	CovGammaAlpha       /= fN*Sgm.fLPowerLawGamma*Sgm.fLPowerLawAlpha;
	CovLogGammaAlpha    /= fN*LogSgm.fLPowerLawGamma*Sgm.fLPowerLawAlpha;
	CovGammaBeta        /= fN*Sgm.fLPowerLawGamma*Sgm.fLPowerLawBeta;
	CovLogGammaBeta     /= fN*LogSgm.fLPowerLawGamma*Sgm.fLPowerLawBeta;
	CovGammaKappa       /= fN*Sgm.fLPowerLawGamma*Sgm.fBDecayKappa;
	CovLogGammaKappa    /= fN*LogSgm.fLPowerLawGamma*Sgm.fBDecayKappa;
	CovGammaDelta       /= fN*Sgm.fLPowerLawGamma*Sgm.fBDecayDelta;
	CovLogGammaDelta    /= fN*LogSgm.fLPowerLawGamma*Sgm.fBDecayDelta;
	CovGammaLogDelta    /= fN*Sgm.fLPowerLawGamma*LogSgm.fBDecayDelta;
	CovLogGammaLogDelta /= fN*LogSgm.fLPowerLawGamma*LogSgm.fBDecayDelta;

	CovAlphaBeta        /= fN*Sgm.fLPowerLawAlpha*Sgm.fLPowerLawBeta;
	CovAlphaKappa       /= fN*Sgm.fLPowerLawAlpha*Sgm.fBDecayKappa;
	CovAlphaDelta       /= fN*Sgm.fLPowerLawAlpha*Sgm.fBDecayDelta;
	CovAlphaLogDelta    /= fN*Sgm.fLPowerLawAlpha*LogSgm.fBDecayDelta;

	CovBetaKappa        /= fN*Sgm.fLPowerLawBeta*Sgm.fBDecayKappa;
	CovBetaDelta        /= fN*Sgm.fLPowerLawBeta*Sgm.fBDecayDelta;
	CovBetaLogDelta     /= fN*Sgm.fLPowerLawBeta*LogSgm.fBDecayDelta;

	CovKappaDelta       /= fN*Sgm.fBDecayKappa*Sgm.fBDecayDelta;
	CovKappaLogDelta    /= fN*Sgm.fBDecayKappa*LogSgm.fBDecayDelta;


	ofstream out2("Results/Avg.dat");
	out2<<"LinGamma "<<Avg.fLPowerLawGamma<<" "<<Sgm.fLPowerLawGamma<<endl;
	out2<<"LinAlpha "<<Avg.fLPowerLawAlpha<<" "<<Sgm.fLPowerLawAlpha<<endl;
	out2<<"LinBeta "<<Avg.fLPowerLawBeta<<" "<<Sgm.fLPowerLawBeta<<endl;
	out2<<"LinDelta "<<Avg.fBDecayDelta<<" "<<Sgm.fBDecayDelta<<endl;
	out2<<"LinKappa "<<Avg.fBDecayKappa<<" "<<Sgm.fBDecayKappa<<endl;
	out2<<endl;
	out2<<"LogGamma "<<LogAvg.fLPowerLawGamma<<" "<<LogSgm.fLPowerLawGamma<<endl;
	out2<<"LogAlpha "<<LogAvg.fLPowerLawAlpha<<" "<<LogSgm.fLPowerLawAlpha<<endl;
	out2<<"LogBeta "<<LogAvg.fLPowerLawBeta<<" "<<LogSgm.fLPowerLawBeta<<endl;
	out2<<"LogDelta "<<LogAvg.fBDecayDelta<<" "<<LogSgm.fBDecayDelta<<endl;
	out2<<"LogKappa "<<LogAvg.fBDecayKappa<<" "<<LogSgm.fBDecayKappa<<endl;
	out2.close();


	ofstream cor("Results/Correlation.tex");

	cor<<"\\begin{table*}"<<endl;
	cor<<"\\begin{tabular}{c  c  c  c  c  c  c  c }"<<endl;
	cor<<" & $\\gamma$ & $\\log_{10}\\gamma$ & $\\alpha$ & $\\beta$ & $\\Delta$ & $\\log_{10}\\Delta$ & $\\kappa$ \\\\"<<endl;
	cor<<"\\hline"<<endl;
	cor<<"\\hline"<<endl;
	cor<<"$\\gamma$ & $1$ & & $"<<CovGammaAlpha<<"$ & $"<< CovGammaBeta<<"$ & $"<<CovGammaDelta<<"$ & $"<<CovGammaLogDelta<<"$ & $"<<CovGammaKappa<<"$\\\\"<<endl;
	cor<<"$\\log_{10}\\gamma$ & & $1$ & $"<<CovLogGammaAlpha<<"$ & $"<< CovLogGammaBeta<<"$ & $"<<CovLogGammaDelta<<"$ & $"<<CovLogGammaLogDelta<<"$ & $"<<CovLogGammaKappa<<"$\\\\"<<endl;
	cor<<"\\hline"<<endl;
	cor<<"$\\alpha$ & & & $1$ & $"<< CovAlphaBeta<<"$ & $"<<CovAlphaDelta<<"$ & $"<<CovAlphaLogDelta<<"$ & $"<<CovAlphaKappa<<"$\\\\"<<endl;
	cor<<"\\hline"<<endl;
	cor<<"$\\beta$ & & & & $1$ & $"<<CovBetaDelta<<"$ & $"<<CovBetaLogDelta<<"$ & $"<<CovBetaKappa<<"$\\\\"<<endl;
	cor<<"\\hline"<<endl;
	cor<<"$\\Delta$ & & & & & $1$ & & $"<<CovKappaDelta<<"$\\\\"<<endl;
	cor<<"$\\log{10}\\Delta$ & & & & & & $1$ & $"<<CovKappaLogDelta<<"$\\\\"<<endl;
	cor<<"\\hline"<<endl;
	cor<<"$\\kappa$ & & & & & & & $1$ \\\\"<<endl;
	cor<<"\\hline"<<endl;
	cor<<"\\hline"<<endl;
	cor<<"\\end{tabular}"<<endl;
	cor<<"\\end{table*}"<<endl;

	cor.close();

#if 0
	if (m_pCfgPtr->ModelConstraints.bBDecayTimeDeltaMean) 
	{
		HistBDecayTimeDeltaMean->Densify();
		HistBDecayTimeDeltaMean->Write("BDecayTimeDeltaMean.hist");
		delete HistBDecayTimeDeltaMean;
	}

	if (m_pCfgPtr->ModelConstraints.bBDecayTimeDeltaSigma) 
	{
		HistBDecayTimeDeltaSigma->Densify();
		HistBDecayTimeDeltaSigma->Write("BDecayTimeDeltaSigma.hist");
		delete HistBDecayTimeDeltaSigma;
	}
#endif
}

void CMCMC::WriteBinary(const char* Filename)
{
	//cout<<"Writing output to: "<<Filename<<endl;
	ofstream out(Filename,std::ofstream::binary);

	if (!out.is_open()) ERROR("Can't open file for writing");


	const size_t nSize = sizeof(SModelParameters);
	//if(nSize != MODEL_PARAMETERS_SIZE) ERROR("Struct size different than predefined");


	//This is not trival so it has to be written
	out.write(reinterpret_cast<char*>(&m_pCfgPtr->Sim.nNumOfChains),sizeof(int));
	out.write(reinterpret_cast<char*>(&m_pCfgPtr->Sim.nNumPerChain),sizeof(int));


	//UModelParameters Buffer;
	for(int i(0); i<m_pCfgPtr->Sim.nNumOfChains; i++)
	{
		for(int j(0); j<m_pCfgPtr->Sim.nNumPerChain; j++)
		{
			//Buffer.Data=m_pChain[i][j];
			out.write(reinterpret_cast<char *>(&(m_pChain[i][j])), nSize);
		}
	}	
	out.close();
}


void CMCMC::WriteBinaryHDF5(const char* Filename) const
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
	double *pfBuf = new double[m_nNumOfChains*m_nNumPerChain];
	int *pnBuf = new int[m_nNumOfChains*m_nNumPerChain];

   // Create HDF5 file and dataset
	H5::H5File HDF5File(Filename, H5F_ACC_TRUNC);
	const int nRank(2);
	hsize_t pnDims[nRank];
	pnDims[1]=static_cast<hsize_t>(m_nNumOfChains);//X
	pnDims[0]=static_cast<hsize_t>(m_nNumPerChain);//Y
    	H5::DataSpace Space(nRank, pnDims);

	H5::DataSet Set;


	Set = HDF5File.createDataSet("nChainNumber", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++) for(int y(0);y<m_nNumPerChain;y++) pnBuf[x + y*m_nNumOfChains] = x; 
	Set.write(pnBuf,H5::PredType::NATIVE_INT);

	Set = HDF5File.createDataSet("nLinkNumber", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++) for(int y(0);y<m_nNumPerChain;y++) pnBuf[x + y*m_nNumOfChains] = y; 
	Set.write(pnBuf,H5::PredType::NATIVE_INT);


	Set = HDF5File.createDataSet("fLPowerLawGamma", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++) for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fLPowerLawGamma; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fLPowerLawAlpha", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fLPowerLawAlpha; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fLPowerLawBeta", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fLPowerLawBeta; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fBDecayDelta", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fBDecayDelta; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fBDecayKappa", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fBDecayKappa; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fBInitMean", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fBInitMean; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fBInitSigma", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fBInitSigma; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fBFinalMean", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fBFinalMean; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fBFinalSigma", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fBFinalSigma; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fPInitMean", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fPInitMean; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fPInitSigma", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fPInitSigma; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fVKickMod1Mean", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fVKickMod1Mean; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fVKickMod1Sigma", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fVKickMod1Sigma; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fVKickMod2Mean", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fVKickMod2Mean; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fVKickMod2Sigma", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fVKickMod2Sigma; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fLnLikelihood", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fLnLikelihood; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fDValue", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fDValue; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	Set = HDF5File.createDataSet("fDeathPsi", H5::PredType::NATIVE_DOUBLE, Space);
	for(int x(0);x<m_nNumOfChains;x++)for(int y(0);y<m_nNumPerChain;y++) pfBuf[x + y*m_nNumOfChains] = m_pChain[x][y].fDeathPsi; 
	Set.write(pfBuf,H5::PredType::NATIVE_DOUBLE);

	delete [] pfBuf;


}


void CMCMC::CombineBin(const char* OutputFilename,  char** InputFilenameVec, const int nNumOfFiles)
{
	//We don't need this 
	m_DeleteChains();

	ofstream out(OutputFilename, std::ofstream::binary);
	int nChains(0), nLinks(0);

//read first
	ifstream in(InputFilenameVec[0],std::ifstream::binary);
	const size_t nSize = sizeof(SModelParameters);

	//if(nSize != MODEL_PARAMETERS_SIZE) ERROR("Struct size different than predefined");

	//get the size of the file
	in.seekg( 0, std::ios::end );
	size_t nFileSize =  in.tellg();
	in.seekg (0, in.beg);

	int nNum=(nFileSize-2*sizeof(int))/nSize;

	if((nFileSize-2*sizeof(int))%nSize > 0) ERROR("Struct size mismatch. Markov chain's file reading not possible.");

	in.read(reinterpret_cast<char*>(&m_pCfgPtr->Sim.nNumOfChains),sizeof(int));
	in.read(reinterpret_cast<char*>(&m_pCfgPtr->Sim.nNumPerChain),sizeof(int));

	if(nNum != (m_pCfgPtr->Sim.nNumOfChains * m_pCfgPtr->Sim.nNumPerChain)) ERROR("Size of the file seem a bit off... Is it corrupted? >:)");

	size_t nToRead = (m_pCfgPtr->Sim.nNumOfChains) * (m_pCfgPtr->Sim.nNumPerChain) * nSize;

	char* pBuffer = new char[nToRead];

	int nTotalChains(m_pCfgPtr->Sim.nNumOfChains * nNumOfFiles);
	out.write(reinterpret_cast<char*>(&nTotalChains),sizeof(int));
	out.write(reinterpret_cast<char*>(&m_pCfgPtr->Sim.nNumPerChain),sizeof(int));


	in.read(pBuffer,nToRead);
	out.write(pBuffer,nToRead);
	in.close();

	for(int i(1);i<nNumOfFiles;i++)
	{
		//assumes same config for each file but checks if true
		in.open(InputFilenameVec[i],std::ifstream::binary);
	
		//get the size of the file
		in.seekg( 0, std::ios::end );
		size_t nFileSize =  in.tellg();
		in.seekg (0, in.beg);

		nNum=(nFileSize-2*sizeof(int))/nSize;

		if((nFileSize-2*sizeof(int))%nSize > 0) ERROR("Struct size mismatch. Markov chain's file reading not possible.");
		in.read(reinterpret_cast<char*>(&nChains),sizeof(int));
		in.read(reinterpret_cast<char*>(&nLinks),sizeof(int));

		if(nNum != (m_pCfgPtr->Sim.nNumOfChains * m_pCfgPtr->Sim.nNumPerChain)) ERROR("Size of the file seem a bit off... Is it corrupted? >:)");
		if( (nChains != m_pCfgPtr->Sim.nNumOfChains) || (nLinks != m_pCfgPtr->Sim.nNumPerChain)) ERROR("Chains geometry different betwean files");

		in.read(pBuffer,nToRead);
		out.write(pBuffer,nToRead);


		in.close();
	}

	out.close();
	delete [] pBuffer;

}




void CMCMC::ReadBinary(const char* Filename)
{
	ifstream in(Filename,std::fstream::binary);
	const size_t nSize = sizeof(SModelParameters);

	//if(nSize != MODEL_PARAMETERS_SIZE) ERROR("Struct size different than predefined");

	m_DeleteChains();

	
	//get the size of the file
	in.seekg( 0, std::ios::end );
	size_t nFileSize =  in.tellg();
	in.seekg (0, in.beg);

	int nNum=(nFileSize-2*sizeof(int))/nSize;
	if((nFileSize-2*sizeof(int))%nSize > 0) ERROR("Struct size mismatch. Markov chain's file reading not possible.");
	in.read(reinterpret_cast<char*>(&m_pCfgPtr->Sim.nNumOfChains),sizeof(int));
	in.read(reinterpret_cast<char*>(&m_pCfgPtr->Sim.nNumPerChain),sizeof(int));
	m_CreateEmptyChains();

	if(nNum != (m_pCfgPtr->Sim.nNumOfChains * m_pCfgPtr->Sim.nNumPerChain)) ERROR("Size of the file seem a bit off... Is it corrupted? >:)");

	//UModelParameters Buffer;
	for(int i(0); i<m_pCfgPtr->Sim.nNumOfChains; i++)
	{
		for(int j(0); j<m_pCfgPtr->Sim.nNumPerChain; j++)
		{
			in.read(reinterpret_cast<char *>(&(m_pChain[i][j])), nSize);
	//		cerr<<i<<" "<<j<<" "<<m_pChain[i][j].fLnLikelihood<<endl;
		}
	}	
	in.close();

	m_nNumOfChains=m_pCfgPtr->Sim.nNumOfChains;
	m_nNumPerChain=m_pCfgPtr->Sim.nNumPerChain;
	//cout<<m_nNumOfChains<<" "<<m_nNumPerChain<<endl;
}







void CMCMCPostProcessing::FindExtremes()
{
	MinValues=MaxValues=m_pChain[0][0];

	for(int i(0);i<m_nNumOfChains;i++)
	{
		for(int j(0);j<m_nNumPerChain;j++)
		{
			if(m_pChain[i][j].fLPowerLawGamma > MaxValues.fLPowerLawGamma)
			{
				MaxValues.fLPowerLawGamma=m_pChain[i][j].fLPowerLawGamma;
			}
			if(m_pChain[i][j].fLPowerLawGamma < MinValues.fLPowerLawGamma)
			{
				MinValues.fLPowerLawGamma=m_pChain[i][j].fLPowerLawGamma;
			}

			if(m_pChain[i][j].fLPowerLawAlpha > MaxValues.fLPowerLawAlpha)
			{
				MaxValues.fLPowerLawAlpha=m_pChain[i][j].fLPowerLawAlpha;
			}
			if(m_pChain[i][j].fLPowerLawAlpha < MinValues.fLPowerLawAlpha)
			{
				MinValues.fLPowerLawAlpha=m_pChain[i][j].fLPowerLawAlpha;
			}

			if(m_pChain[i][j].fLPowerLawBeta > MaxValues.fLPowerLawBeta)
			{
				MaxValues.fLPowerLawBeta=m_pChain[i][j].fLPowerLawBeta;
			}
			if(m_pChain[i][j].fLPowerLawBeta < MinValues.fLPowerLawBeta)
			{
				MinValues.fLPowerLawBeta=m_pChain[i][j].fLPowerLawBeta;
			}

			if(m_pChain[i][j].fBDecayDelta > MaxValues.fBDecayDelta)
			{
				MaxValues.fBDecayDelta=m_pChain[i][j].fBDecayDelta;
			}
			if(m_pChain[i][j].fBDecayDelta < MinValues.fBDecayDelta)
			{
				MinValues.fBDecayDelta=m_pChain[i][j].fBDecayDelta;
			}

			if(m_pChain[i][j].fBDecayKappa > MaxValues.fBDecayKappa)
			{
				MaxValues.fBDecayKappa=m_pChain[i][j].fBDecayKappa;
			}
			if(m_pChain[i][j].fBDecayKappa < MinValues.fBDecayKappa)
			{
				MinValues.fBDecayKappa=m_pChain[i][j].fBDecayKappa;
			}

	

			if(m_pChain[i][j].fLnLikelihood > MaxValues.fLnLikelihood)
			{
				MaxValues.fLnLikelihood=m_pChain[i][j].fLnLikelihood;
			}
			if(m_pChain[i][j].fLnLikelihood < MinValues.fLnLikelihood)
			{
				MinValues.fLnLikelihood=m_pChain[i][j].fLnLikelihood;
			}
		}
	}

	MaxValues.fLPowerLawGamma=static_cast<int>(log10(MaxValues.fLPowerLawGamma))+1.;
	MaxValues.fLPowerLawAlpha=static_cast<int>(MaxValues.fLPowerLawAlpha)+1.;
	MaxValues.fLPowerLawBeta=static_cast<int>(MaxValues.fLPowerLawBeta)+1.;
	MaxValues.fBDecayDelta=static_cast<int>(log10(MaxValues.fBDecayDelta))+1.;
	MaxValues.fBDecayKappa=static_cast<int>(MaxValues.fBDecayKappa)+1.;

	MinValues.fLPowerLawGamma=static_cast<int>(log10(MinValues.fLPowerLawGamma))-1.;
	MinValues.fLPowerLawAlpha=static_cast<int>(MinValues.fLPowerLawAlpha)-1.;
	MinValues.fLPowerLawBeta=static_cast<int>(MinValues.fLPowerLawBeta)-1.;
//TODO Delta ma bledna wartosc
	MinValues.fBDecayDelta=static_cast<int>(log10(MinValues.fBDecayDelta))-1.;
	MinValues.fBDecayKappa=static_cast<int>(MinValues.fBDecayKappa)-1.;



}
void CMCMCPostProcessing::WriteExtremes()
{
	ofstream out;
	out.open("Results/Min.txt");

	out<<MinValues.fLPowerLawGamma<<" "
	   <<MinValues.fLPowerLawAlpha<<" "
	   <<MinValues.fLPowerLawBeta<<" "
	   <<MinValues.fBDecayDelta<<" "
	   <<MinValues.fBDecayKappa<<" "
	   <<MinValues.fLnLikelihood<<endl;

	out.close();


	out.open("Results/Max.txt");
	out<<MaxValues.fLPowerLawGamma<<" "
	   <<MaxValues.fLPowerLawAlpha<<" "
	   <<MaxValues.fLPowerLawBeta<<" "
	   <<MaxValues.fBDecayDelta<<" "
	   <<MaxValues.fBDecayKappa<<" "
	   <<MaxValues.fLnLikelihood<<endl;
	out.close();


}

void CMCMCPostProcessing::SetExtremes()
{
	MaxValues.fLPowerLawGamma=static_cast<int>(log10(m_pCfgPtr->ModelConstraints.Max.fLPowerLawGamma))+1;
	MaxValues.fLPowerLawAlpha=static_cast<int>((m_pCfgPtr->ModelConstraints.Max.fLPowerLawAlpha))+1;
	MaxValues.fLPowerLawBeta=static_cast<int>((m_pCfgPtr->ModelConstraints.Max.fLPowerLawBeta))+1;
	MaxValues.fBDecayDelta=static_cast<int>(log10(m_pCfgPtr->ModelConstraints.Max.fBDecayDelta))+1.;
	MaxValues.fBDecayKappa=static_cast<int>( (m_pCfgPtr->ModelConstraints.Max.fBDecayKappa))+1;

	MinValues.fLPowerLawGamma=static_cast<int>(log10(m_pCfgPtr->ModelConstraints.Min.fLPowerLawGamma))-1;
	MinValues.fLPowerLawAlpha=static_cast<int>((m_pCfgPtr->ModelConstraints.Min.fLPowerLawAlpha))-1;
	MinValues.fLPowerLawBeta=static_cast<int>((m_pCfgPtr->ModelConstraints.Min.fLPowerLawBeta))-1;
	MinValues.fBDecayDelta=static_cast<int>(log10(m_pCfgPtr->ModelConstraints.Min.fBDecayDelta))-1.;
	MinValues.fBDecayKappa=static_cast<int>( (m_pCfgPtr->ModelConstraints.Min.fBDecayKappa))-1;

}

void CMCMCPostProcessing::MakeHistograms(const int nSquareResolution)
{
	FindExtremes();
//	SetExtremes();
	WriteExtremes();
	stringstream translator;

	ofstream out("Results/Res.txt");
	out<<nSquareResolution<<endl;
	out.close();

	CAreaDensity AlphaBeta(MinValues.fLPowerLawAlpha, MinValues.fLPowerLawBeta, MaxValues.fLPowerLawAlpha, MaxValues.fLPowerLawBeta, nSquareResolution, nSquareResolution);
	CAreaDensity AlphaGamma(MinValues.fLPowerLawAlpha, MinValues.fLPowerLawGamma, MaxValues.fLPowerLawAlpha, MaxValues.fLPowerLawGamma, nSquareResolution, nSquareResolution);

	CAreaDensity AlphaDelta(MinValues.fLPowerLawAlpha, MinValues.fBDecayDelta, MaxValues.fLPowerLawAlpha, MaxValues.fBDecayDelta, nSquareResolution, nSquareResolution);
	CAreaDensity AlphaKappa(MinValues.fLPowerLawAlpha, MinValues.fBDecayKappa, MaxValues.fLPowerLawAlpha, MaxValues.fBDecayKappa, nSquareResolution, nSquareResolution);
	
	CAreaDensity BetaGamma(MinValues.fLPowerLawBeta, MinValues.fLPowerLawGamma, MaxValues.fLPowerLawBeta, MaxValues.fLPowerLawGamma, nSquareResolution, nSquareResolution);
	CAreaDensity BetaDelta(MinValues.fLPowerLawBeta, MinValues.fBDecayDelta, MaxValues.fLPowerLawBeta, MaxValues.fBDecayDelta, nSquareResolution, nSquareResolution);
	CAreaDensity BetaKappa(MinValues.fLPowerLawBeta, MinValues.fBDecayKappa, MaxValues.fLPowerLawBeta, MaxValues.fBDecayKappa, nSquareResolution, nSquareResolution);
	CAreaDensity DeltaKappa(MinValues.fBDecayDelta, MinValues.fBDecayKappa, MaxValues.fBDecayDelta, MaxValues.fBDecayKappa, nSquareResolution, nSquareResolution);
	CAreaDensity DeltaGamma(MinValues.fBDecayDelta, MinValues.fLPowerLawGamma, MaxValues.fBDecayDelta, MaxValues.fLPowerLawGamma, nSquareResolution, nSquareResolution);
	CAreaDensity KappaGamma(MinValues.fBDecayKappa, MinValues.fLPowerLawGamma, MaxValues.fBDecayKappa, MaxValues.fLPowerLawGamma, nSquareResolution, nSquareResolution);

	int nStep=5;//przesunąć do configu co by się mogło stać parametrem i
	//for(int nPointOfChain(0);nPointOfChain<1;nPointOfChain+=nStep)
	for(int nPointOfChain(0);nPointOfChain<m_nNumPerChain;nPointOfChain+=nStep)
	{
		//cout<<nPointOfChain<<" of "<<m_nNumPerChain<<" every "<<nStep<<endl;
		for(int nChain(0);nChain<m_nNumOfChains;nChain++)
		{
			AlphaBeta.AddPoint(m_pChain[nChain][nPointOfChain].fLPowerLawAlpha, m_pChain[nChain][nPointOfChain].fLPowerLawBeta);
			AlphaGamma.AddPoint(m_pChain[nChain][nPointOfChain].fLPowerLawAlpha,log10( m_pChain[nChain][nPointOfChain].fLPowerLawGamma));
			AlphaDelta.AddPoint(m_pChain[nChain][nPointOfChain].fLPowerLawAlpha,log10( m_pChain[nChain][nPointOfChain].fBDecayDelta));
			AlphaKappa.AddPoint(m_pChain[nChain][nPointOfChain].fLPowerLawAlpha, m_pChain[nChain][nPointOfChain].fBDecayKappa);
			BetaGamma.AddPoint(m_pChain[nChain][nPointOfChain].fLPowerLawBeta,log10( m_pChain[nChain][nPointOfChain].fLPowerLawGamma));
			BetaDelta.AddPoint(m_pChain[nChain][nPointOfChain].fLPowerLawBeta,log10( m_pChain[nChain][nPointOfChain].fBDecayDelta));
			BetaKappa.AddPoint(m_pChain[nChain][nPointOfChain].fLPowerLawBeta, m_pChain[nChain][nPointOfChain].fBDecayKappa);
			DeltaKappa.AddPoint(log10(m_pChain[nChain][nPointOfChain].fBDecayDelta), m_pChain[nChain][nPointOfChain].fBDecayKappa);
			DeltaGamma.AddPoint(log10(m_pChain[nChain][nPointOfChain].fBDecayDelta), log10(m_pChain[nChain][nPointOfChain].fLPowerLawGamma));
			KappaGamma.AddPoint(m_pChain[nChain][nPointOfChain].fBDecayKappa,log10( m_pChain[nChain][nPointOfChain].fLPowerLawGamma));
			
			//cout<<log10(m_pChain[nChain][nPointOfChain].fBDecayDelta)<<" "<< m_pChain[nChain][nPointOfChain].fBDecayKappa<<endl;
		}

		AlphaBeta.Norm();
		AlphaGamma.Norm();
		AlphaDelta.Norm();
		AlphaKappa.Norm();
		BetaGamma.Norm();
		BetaDelta.Norm();
		BetaKappa.Norm();
		DeltaKappa.Norm();
		DeltaGamma.Norm();
		KappaGamma.Norm();

		AlphaBeta.Sigmify();
		AlphaGamma.Sigmify();
		AlphaDelta.Sigmify();
		AlphaKappa.Sigmify();
		BetaGamma.Sigmify();
		BetaDelta.Sigmify();
		BetaKappa.Sigmify();
		DeltaKappa.Sigmify();
		DeltaGamma.Sigmify();
		KappaGamma.Sigmify();

		translator.str("");
		translator.clear();
		translator<<"Results/AlphaBeta."<<nPointOfChain<<".mtx";
		AlphaBeta.Write(translator.str().c_str());
		
		translator.str("");
		translator.clear();
		translator<<"Results/AlphaGamma."<<nPointOfChain<<".mtx";
		AlphaGamma.Write(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/AlphaDelta."<<nPointOfChain<<".mtx";
		AlphaDelta.Write(translator.str().c_str());


		translator.str("");
		translator.clear();
		translator<<"Results/AlphaKappa."<<nPointOfChain<<".mtx";
		AlphaKappa.Write(translator.str().c_str());
		
		translator.str("");
		translator.clear();
		translator<<"Results/BetaGamma."<<nPointOfChain<<".mtx";
		BetaGamma.Write(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/BetaDelta."<<nPointOfChain<<".mtx";
		BetaDelta.Write(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/BetaKappa."<<nPointOfChain<<".mtx";
		BetaKappa.Write(translator.str().c_str());


		translator.str("");
		translator.clear();
		translator<<"Results/DeltaKappa."<<nPointOfChain<<".mtx";
		DeltaKappa.Write(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/DeltaGamma."<<nPointOfChain<<".mtx";
		DeltaGamma.Write(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/KappaGamma."<<nPointOfChain<<".mtx";
		KappaGamma.Write(translator.str().c_str());






		translator.str("");
		translator.clear();
		translator<<"Results/AlphaBeta."<<nPointOfChain<<".sigma.mtx";
		AlphaBeta.WriteSigma(translator.str().c_str());
		
		translator.str("");
		translator.clear();
		translator<<"Results/AlphaGamma."<<nPointOfChain<<".sigma.mtx";
		AlphaGamma.WriteSigma(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/AlphaDelta."<<nPointOfChain<<".sigma.mtx";
		AlphaDelta.WriteSigma(translator.str().c_str());


		translator.str("");
		translator.clear();
		translator<<"Results/AlphaKappa."<<nPointOfChain<<".sigma.mtx";
		AlphaKappa.WriteSigma(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/BetaGamma."<<nPointOfChain<<".sigma.mtx";
		BetaGamma.WriteSigma(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/BetaDelta."<<nPointOfChain<<".sigma.mtx";
		BetaDelta.WriteSigma(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/BetaKappa."<<nPointOfChain<<".sigma.mtx";
		BetaKappa.WriteSigma(translator.str().c_str());


		translator.str("");
		translator.clear();
		translator<<"Results/DeltaKappa."<<nPointOfChain<<".sigma.mtx";
		DeltaKappa.WriteSigma(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/DeltaGamma."<<nPointOfChain<<".sigma.mtx";
		DeltaGamma.WriteSigma(translator.str().c_str());

		translator.str("");
		translator.clear();
		translator<<"Results/KappaGamma."<<nPointOfChain<<".sigma.mtx";
		KappaGamma.WriteSigma(translator.str().c_str());



		AlphaBeta.Reset();
		AlphaGamma.Reset();
		AlphaDelta.Reset();
		AlphaKappa.Reset();
		BetaGamma.Reset();
		BetaDelta.Reset();
		BetaKappa.Reset();
		DeltaKappa.Reset();
		DeltaGamma.Reset();
		KappaGamma.Reset();

	}
}






///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

void CModelParameterDistributions::ScaleJumpSigmas()
{
    if(m_pCfgPtr->Sim.bMCMCVerbose) cout<<m_pCfgPtr->Sim.fJumpScallingFactor<<endl;
	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawGamma)
	{
		m_pCfgPtr->MCMCStepSize.fLPowerLawGamma = pow(10., log10(m_pCfgPtr->MCMCStepSize.fLPowerLawGamma)* m_pCfgPtr->Sim.fJumpScallingFactor);
	}

	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawAlpha)
	{
		m_pCfgPtr->MCMCStepSize.fLPowerLawAlpha *=m_pCfgPtr->Sim.fJumpScallingFactor;
	}

	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawBeta)
	{
		m_pCfgPtr->MCMCStepSize.fLPowerLawBeta *=m_pCfgPtr->Sim.fJumpScallingFactor;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBDecayDelta)
	{
		m_pCfgPtr->MCMCStepSize.fBDecayDelta  = pow(10., log10(m_pCfgPtr->MCMCStepSize.fBDecayDelta)* 0.5);
	}

	if(m_pCfgPtr->MCMCParameterUse.bBDecayKappa)
	{
		m_pCfgPtr->MCMCStepSize.fBDecayKappa *=m_pCfgPtr->Sim.fJumpScallingFactor;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBInitMean)
	{
		m_pCfgPtr->MCMCStepSize.fBInitMean *=m_pCfgPtr->Sim.fJumpScallingFactor;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBInitSigma)
	{
		m_pCfgPtr->MCMCStepSize.fBInitSigma *=m_pCfgPtr->Sim.fJumpScallingFactor;
	}

	if(m_pCfgPtr->MCMCParameterUse.bPInitMean)
	{
		m_pCfgPtr->MCMCStepSize.fPInitMean *=m_pCfgPtr->Sim.fJumpScallingFactor;
	}

	if(m_pCfgPtr->MCMCParameterUse.bPInitSigma)
	{
		m_pCfgPtr->MCMCStepSize.fPInitSigma *=m_pCfgPtr->Sim.fJumpScallingFactor;
	}

	if(m_pCfgPtr->MCMCParameterUse.bDeathPsi)
	{
		m_pCfgPtr->MCMCStepSize.fDeathPsi *=m_pCfgPtr->Sim.fJumpScallingFactor;
	}




}






SModelParameters CModelParameterDistributions::InitialDistribution() const
{
	SModelParameters Result;

	//Start by filling the Result with default parameters
	Result = m_pCfgPtr->Model;

	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawGamma)
	{
		Result.fLPowerLawGamma = m_PowerLawGammaInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fLPowerLawGamma: "<<Result.fLPowerLawGamma<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawAlpha)
	{
		Result.fLPowerLawAlpha = m_PowerLawAlphaInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fLPowerLawAlpha: "<<Result.fLPowerLawAlpha<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawBeta)
	{
		Result.fLPowerLawBeta = m_PowerLawBetaInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fLPowerLawBeta: "<<Result.fLPowerLawBeta<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBDecayDelta)
	{
		Result.fBDecayDelta = m_BDecayDeltaInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBDecayDelta: "<<Result.fBDecayDelta<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBDecayKappa)
	{
		Result.fBDecayKappa = m_BDecayKappaInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBDecayKappa: "<<Result.fBDecayKappa<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBInitMean)
	{
		Result.fBInitMean = m_BInitMeanInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBInitMean: "<<Result.fBInitMean<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBInitSigma)
	{
		Result.fBInitSigma = m_BInitSigmaInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBInitSigma: "<<Result.fBInitSigma<<endl;
	}



	if(m_pCfgPtr->MCMCParameterUse.bPInitMean)
	{
		Result.fPInitMean = m_PInitMeanInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fPInitMean: "<<Result.fPInitMean<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bPInitSigma)
	{
		Result.fPInitSigma = m_PInitSigmaInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fPInitSigma: "<<Result.fPInitSigma<<endl;
	}


	if(m_pCfgPtr->MCMCParameterUse.bDeathPsi)
	{
		Result.fDeathPsi = m_DeathPsiInitialDistribution(); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fDeathPsi: "<<Result.fDeathPsi<<endl;
	}


	return Result;
}

SModelParameters CModelParameterDistributions::SamplingDistribution(const SModelParameters Current) const
{
	SModelParameters Result;
	Result = m_pCfgPtr->Model;
    Result.fLnLikelihood=0.;

	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawGamma)
	{
		Result.fLPowerLawGamma = m_PowerLawGammaSamplingDistribution(Current.fLPowerLawGamma); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fLPowerLawGamma: "<<Result.fLPowerLawGamma<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawAlpha)
	{
		Result.fLPowerLawAlpha = m_PowerLawAlphaSamplingDistribution(Current.fLPowerLawAlpha); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fLPowerLawAlpha: "<<Result.fLPowerLawAlpha<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bLPowerLawBeta)
	{
		Result.fLPowerLawBeta = m_PowerLawBetaSamplingDistribution(Current.fLPowerLawBeta); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fLPowerLawBeta: "<<Result.fLPowerLawBeta<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBDecayDelta)
	{
		Result.fBDecayDelta = m_BDecayDeltaSamplingDistribution(Current.fBDecayDelta); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBDecayDelta: "<<Result.fBDecayDelta<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBDecayKappa)
	{
		Result.fBDecayKappa = m_BDecayKappaSamplingDistribution(Current.fBDecayKappa); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBDecayKappa: "<<Result.fBDecayKappa<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBDecayKappa)
	{
		Result.fBDecayKappa = m_BDecayKappaSamplingDistribution(Current.fBDecayKappa); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBDecayKappa: "<<Result.fBDecayKappa<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBInitMean)
	{
		Result.fBInitMean = m_BInitMeanSamplingDistribution(Current.fBInitMean); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBInitMean: "<<Result.fBInitMean<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bBInitSigma)
	{
		Result.fBInitSigma = m_BInitSigmaSamplingDistribution(Current.fBInitSigma); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fBInitSigma: "<<Result.fBInitSigma<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bPInitMean)
	{
		Result.fPInitMean = m_PInitMeanSamplingDistribution(Current.fPInitMean); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fPInitMean: "<<Result.fPInitMean<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bPInitSigma)
	{
		Result.fPInitSigma = m_PInitSigmaSamplingDistribution(Current.fPInitSigma); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fPInitSigma: "<<Result.fPInitSigma<<endl;
	}

	if(m_pCfgPtr->MCMCParameterUse.bDeathPsi)
	{
		Result.fDeathPsi = m_DeathPsiSamplingDistribution(Current.fDeathPsi); 
		if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) cout<<"Result.fDeathPsi: "<<Result.fDeathPsi<<endl;
	}

	return Result;
}


RealType CModelParameterDistributions::m_PowerLawGammaInitialDistribution() const
{
	return pow(10.,m_pRanPtr->FlatMinMax( log10(m_pCfgPtr->ModelConstraints.Min.fLPowerLawGamma), log10(m_pCfgPtr->ModelConstraints.Max.fLPowerLawGamma) ) ) ;	
} 

RealType CModelParameterDistributions::m_PowerLawAlphaInitialDistribution() const
{
	return m_pRanPtr->FlatMinMax(m_pCfgPtr->ModelConstraints.Min.fLPowerLawAlpha, m_pCfgPtr->ModelConstraints.Max.fLPowerLawAlpha);
}

RealType CModelParameterDistributions::m_PowerLawBetaInitialDistribution()  const
{
	return m_pRanPtr->FlatMinMax(m_pCfgPtr->ModelConstraints.Min.fLPowerLawBeta, m_pCfgPtr->ModelConstraints.Max.fLPowerLawBeta);
}

//BDecay
RealType CModelParameterDistributions::m_BDecayDeltaInitialDistribution() const
{
	return pow(10., m_pRanPtr->FlatMinMax( log10(m_pCfgPtr->ModelConstraints.Min.fBDecayDelta), log10(m_pCfgPtr->ModelConstraints.Max.fBDecayDelta) ) );
}

RealType CModelParameterDistributions::m_BDecayKappaInitialDistribution() const
{
	return m_pRanPtr->FlatMinMax(m_pCfgPtr->ModelConstraints.Min.fBDecayKappa, m_pCfgPtr->ModelConstraints.Max.fBDecayKappa);
}

RealType CModelParameterDistributions::m_BInitMeanInitialDistribution() const
{
	return m_pRanPtr->FlatMinMax(m_pCfgPtr->ModelConstraints.Min.fBInitMean, m_pCfgPtr->ModelConstraints.Max.fBInitMean);
}

RealType CModelParameterDistributions::m_BInitSigmaInitialDistribution() const
{
	return m_pRanPtr->FlatMinMax(m_pCfgPtr->ModelConstraints.Min.fBInitSigma, m_pCfgPtr->ModelConstraints.Max.fBInitSigma);
}

RealType CModelParameterDistributions::m_PInitMeanInitialDistribution() const
{
	return m_pRanPtr->FlatMinMax(m_pCfgPtr->ModelConstraints.Min.fPInitMean, m_pCfgPtr->ModelConstraints.Max.fPInitMean);
}

RealType CModelParameterDistributions::m_PInitSigmaInitialDistribution() const
{
	return m_pRanPtr->FlatMinMax(m_pCfgPtr->ModelConstraints.Min.fPInitSigma, m_pCfgPtr->ModelConstraints.Max.fPInitSigma);
}

RealType CModelParameterDistributions::m_DeathPsiInitialDistribution() const
{
	return m_pRanPtr->FlatMinMax(m_pCfgPtr->ModelConstraints.Min.fDeathPsi, m_pCfgPtr->ModelConstraints.Max.fDeathPsi);
}

//If you started thinking "lisp macros" ...
//Well yeah, I couldn't agree more. ;]

//PowerLaw
RealType CModelParameterDistributions::m_PowerLawGammaSamplingDistribution(const RealType CurrentGamma) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"PowerLawGamma SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bGammaFail(true);
	int nDraws(0);
	do
	{
		Result = pow( 10., m_pRanPtr->Gauss( log10(m_pCfgPtr->MCMCStepSize.fLPowerLawGamma), log10(CurrentGamma)  ) );

		if(m_pCfgPtr->MCMCParameterConstraint.bLPowerLawGamma)
		{
			bGammaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fLPowerLawGamma) || (Result > m_pCfgPtr->ModelConstraints.Max.fLPowerLawGamma);
		}
		else bGammaFail = false;
		nDraws++;
	}
	while(bGammaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fLPowerLawGamma<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fLPowerLawGamma <<" draws: "<<nDraws<<endl;
	}
	return Result;
}

RealType CModelParameterDistributions::m_PowerLawAlphaSamplingDistribution(const RealType CurrentAlpha) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"PowerLawAlpha SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bAlphaFail(true);
	int nDraws(0);
	do
	{
		Result =   m_pRanPtr->Gauss(m_pCfgPtr->MCMCStepSize.fLPowerLawAlpha, CurrentAlpha  );

		if(m_pCfgPtr->MCMCParameterConstraint.bLPowerLawAlpha)
		{
			bAlphaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fLPowerLawAlpha) || (Result > m_pCfgPtr->ModelConstraints.Max.fLPowerLawAlpha);
		}
		else bAlphaFail = false;
		nDraws++;
	}
	while(bAlphaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fLPowerLawAlpha<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fLPowerLawAlpha <<" draws: "<<nDraws<<endl;
	}

	return Result;
}

RealType CModelParameterDistributions::m_PowerLawBetaSamplingDistribution (const RealType CurrentBeta)  const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"PowerLawBeta SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bBetaFail(true);
	int nDraws(0);
	do
	{
		Result =   m_pRanPtr->Gauss(m_pCfgPtr->MCMCStepSize.fLPowerLawBeta, CurrentBeta );

		if(m_pCfgPtr->MCMCParameterConstraint.bLPowerLawBeta)
		{
			bBetaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fLPowerLawBeta) || (Result > m_pCfgPtr->ModelConstraints.Max.fLPowerLawBeta);
		}
		else bBetaFail = false;
		nDraws++;
	}
	while(bBetaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fLPowerLawBeta<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fLPowerLawBeta <<" draws: "<<nDraws<<endl;
	}
	return Result;
}

//BDecay
//In fact Julia has pretty good macros too.
RealType CModelParameterDistributions::m_BDecayDeltaSamplingDistribution(const RealType CurrentDelta) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"BDecayDelta SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bDeltaFail(true);
	int nDraws(0);
	do
	{

		Result =  pow(10.,  m_pRanPtr->Gauss(log10(m_pCfgPtr->MCMCStepSize.fBDecayDelta), log10(CurrentDelta)  ) );

		if(m_pCfgPtr->MCMCParameterConstraint.bBDecayDelta)
		{
			bDeltaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fBDecayDelta) || (Result > m_pCfgPtr->ModelConstraints.Max.fBDecayDelta);
		}
		else bDeltaFail = false;
		nDraws++;

	}
	while(bDeltaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fBDecayDelta<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fBDecayDelta <<" draws: "<<nDraws<<endl;
	}

	return Result;
}


//Damn... I even used regexp's to write those down.
//Such a waste of macro-potential... 
RealType CModelParameterDistributions::m_BDecayKappaSamplingDistribution(const RealType CurrentKappa) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"BDecayKappa SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bKappaFail(true);
	int nDraws(0);
	do
	{
		Result = m_pRanPtr->Gauss(m_pCfgPtr->MCMCStepSize.fBDecayKappa, CurrentKappa );

		if(m_pCfgPtr->MCMCParameterConstraint.bBDecayKappa)
		{
			bKappaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fBDecayKappa) || (Result > m_pCfgPtr->ModelConstraints.Max.fBDecayKappa);
		}
		else bKappaFail = false;
		nDraws++;

	}
	while(bKappaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fBDecayKappa<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fBDecayKappa <<" draws: "<<nDraws<<endl;
	}

	return Result;
}

RealType CModelParameterDistributions::m_BInitMeanSamplingDistribution(const RealType CurrentMean) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"BInitMean SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bKappaFail(true);
	int nDraws(0);
	do
	{
		Result = m_pRanPtr->Gauss(m_pCfgPtr->MCMCStepSize.fBInitMean, CurrentMean );

		if(m_pCfgPtr->MCMCParameterConstraint.bBInitMean)
		{
			bKappaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fBInitMean) || (Result > m_pCfgPtr->ModelConstraints.Max.fBInitMean);
		}
		else bKappaFail = false;
		nDraws++;

	}
	while(bKappaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fBInitMean<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fBInitMean <<" draws: "<<nDraws<<endl;
	}

	return Result;
}

RealType CModelParameterDistributions::m_BInitSigmaSamplingDistribution(const RealType CurrentSigma) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"BInitSigma SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bKappaFail(true);
	int nDraws(0);
	do
	{
		Result = m_pRanPtr->Gauss(m_pCfgPtr->MCMCStepSize.fBInitSigma, CurrentSigma );

		if(m_pCfgPtr->MCMCParameterConstraint.bBInitSigma)
		{
			bKappaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fBInitSigma) || (Result > m_pCfgPtr->ModelConstraints.Max.fBInitSigma);
		}
		else bKappaFail = false;
		nDraws++;

	}
	while(bKappaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fBInitSigma<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fBInitSigma <<" draws: "<<nDraws<<endl;
	}

	return Result;
}

RealType CModelParameterDistributions::m_PInitMeanSamplingDistribution(const RealType CurrentMean) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"PInitMean SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bKappaFail(true);
	int nDraws(0);
	do
	{
		Result = m_pRanPtr->Gauss(m_pCfgPtr->MCMCStepSize.fPInitMean, CurrentMean );

		if(m_pCfgPtr->MCMCParameterConstraint.bPInitMean)
		{
			bKappaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fPInitMean) || (Result > m_pCfgPtr->ModelConstraints.Max.fPInitMean);
		}
		else bKappaFail = false;
		nDraws++;

	}
	while(bKappaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fPInitMean<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fPInitMean <<" draws: "<<nDraws<<endl;
	}

	return Result;
}

RealType CModelParameterDistributions::m_PInitSigmaSamplingDistribution(const RealType CurrentSigma) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"PInitSigma SamplingDistribution"<<endl;
	}
	RealType Result;
	bool bKappaFail(true);
	int nDraws(0);
	do
	{
		Result = m_pRanPtr->Gauss(m_pCfgPtr->MCMCStepSize.fPInitSigma, CurrentSigma );

		if(m_pCfgPtr->MCMCParameterConstraint.bPInitSigma)
		{
			bKappaFail = (Result < m_pCfgPtr->ModelConstraints.Min.fPInitSigma) || (Result > m_pCfgPtr->ModelConstraints.Max.fPInitSigma);
		}
		else bKappaFail = false;
		nDraws++;
		//cout<<nDraws<<" "<<CurrentSigma<<" "<<m_pCfgPtr->MCMCStepSize.fPInitSigma<<" "<<Result<<" "<<bKappaFail<<" "<< m_pCfgPtr->ModelConstraints.Min.fPInitSigma<<" "<<m_pCfgPtr->ModelConstraints.Max.fPInitSigma<<endl;
	}
	while(bKappaFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fPInitSigma<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fPInitSigma <<" draws: "<<nDraws<<endl;
	}

	return Result;
}

RealType CModelParameterDistributions::m_DeathPsiSamplingDistribution(const RealType CurrentPsi) const
{
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"DeathPsi SamplingDistribution"<<endl;
	}

	RealType Result;
	bool bFail(true);
	int nDraws(0);
	do
	{
		Result = m_pRanPtr->Gauss(m_pCfgPtr->MCMCStepSize.fDeathPsi, CurrentPsi );

		if(m_pCfgPtr->MCMCParameterConstraint.bDeathPsi)
		{
			bFail = (Result < m_pCfgPtr->ModelConstraints.Min.fDeathPsi) || (Result > m_pCfgPtr->ModelConstraints.Max.fDeathPsi);
		}
		else bFail = false;
		nDraws++;
	}
	while(bFail);
	if(m_pCfgPtr->Sim.bVerboseParameterUsage && m_pCfgPtr->Sim.bMCMCVerbose) 
	{
		cout<<"min: "<< m_pCfgPtr->ModelConstraints.Min.fDeathPsi<<" max: "<< m_pCfgPtr->ModelConstraints.Max.fDeathPsi <<" draws: "<<nDraws<<endl;
	}

	return Result;
}

void CMCMC::m_ComputeValue(SModelParameters& Parameters)
{
	m_pCfgPtr->Model = Parameters;///< set the parameters for computation

	m_pModel->FullCycle();

	Parameters.fLnLikelihood = m_pCfgPtr->Model.fLnLikelihood;
	Parameters.fDValue = m_pCfgPtr->Model.fDValue;
    if(m_pCfgPtr->Sim.bMCMCVerbose)
    {
    	cout<<"Computed values:"<<endl;	
    	cout<<"\tfLnLikelihood = "<<Parameters.fLnLikelihood<<endl;
    	cout<<"\tfDValue       = "<<Parameters.fDValue<<endl;
    }
}


