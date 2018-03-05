#include<ModelParameters.hpp>
#include<Pulsar.hpp>
#include<cmath>


std::string DescribeModelParameters()
{
	return std::string("#log10(fBDecayDelta) fBDecayKappa log10(fLPowerLawGamma) fLPowerLawAlpha fLPowerLawBeta (...) fLnLikelihood");


}

std::ostream& operator<<(std::ostream& out, const SModelParameters& Model)
{
	out<<log10(Model.fBDecayDelta)
	   <<" "<<Model.fBDecayKappa
	   <<" "<<log10(Model.fLPowerLawGamma)
	   <<" "<<Model.fLPowerLawAlpha
	   <<" "<<Model.fLPowerLawBeta
	   <<" "<<Model.fBInitMean	
	   <<" "<<Model.fBInitSigma	
	   <<" "<<Model.fPInitMean
	   <<" "<<Model.fPInitSigma
	   <<" "<<Model.fVKickMod1Mean
	   <<" "<<Model.fVKickMod1Sigma
	   <<" "<<Model.fVKickMod2Mean
	   <<" "<<Model.fVKickMod2Sigma
	   <<" "<<Model.fDeathPsi
	   <<" "<<Model.fLnLikelihood
	   <<" "<<Model.fDValue
	   <<" "<<Model.nSeed;


	return out;
}

SModelParameters pow(const SModelParameters A, const RealType fB)
{
	SModelParameters Result;
	Result.fBDecayDelta = pow((A.fBDecayDelta), fB);
	Result.fBDecayKappa = pow(A.fBDecayKappa, fB);
	Result.fLPowerLawAlpha = pow(A.fLPowerLawAlpha, fB);
	Result.fLPowerLawBeta = pow(A.fLPowerLawBeta, fB);
	Result.fLPowerLawGamma = pow(A.fLPowerLawGamma, fB);

	Result.fBInitMean = pow(A.fBInitMean, fB);	
	Result.fBInitSigma = pow(A.fBInitSigma, fB);	
	Result.fPInitMean = pow(A.fPInitMean, fB);
	Result.fPInitSigma = pow(A.fPInitSigma, fB);
	Result.fVKickMod1Mean = pow(A.fVKickMod1Mean, fB);
	Result.fVKickMod1Sigma = pow(A.fVKickMod1Sigma, fB);
	Result.fVKickMod2Mean = pow(A.fVKickMod2Mean, fB);
	Result.fVKickMod2Sigma = pow(A.fVKickMod2Sigma, fB);
	Result.fDeathPsi = pow(A.fDeathPsi, fB);


	return Result;	
}

SModelParameters log10(const SModelParameters A)
{
	SModelParameters Result;
	Result.fBDecayDelta = log10(A.fBDecayDelta);
	Result.fBDecayKappa = log10(A.fBDecayKappa);
	Result.fLPowerLawAlpha = log10(A.fLPowerLawAlpha);
	Result.fLPowerLawBeta = log10(A.fLPowerLawBeta);
	Result.fLPowerLawGamma = log10(A.fLPowerLawGamma);
	Result.fBInitMean = log10(A.fBInitMean);	
	Result.fBInitSigma = log10(A.fBInitSigma);	
	Result.fPInitMean = log10(A.fPInitMean);
	Result.fPInitSigma = log10(A.fPInitSigma);
	Result.fVKickMod1Mean = log10(A.fVKickMod1Mean);
	Result.fVKickMod1Sigma = log10(A.fVKickMod1Sigma);
	Result.fVKickMod2Mean = log10(A.fVKickMod2Mean);
	Result.fVKickMod2Sigma = log10(A.fVKickMod2Sigma);
	Result.fDeathPsi = log10(A.fDeathPsi);


	return Result;	

}


SModelParameters sqrt(const SModelParameters A)
{
	SModelParameters Result;
	Result.fBDecayDelta = sqrt(A.fBDecayDelta);
	Result.fBDecayKappa = sqrt(A.fBDecayKappa);
	Result.fLPowerLawAlpha = sqrt(A.fLPowerLawAlpha);
	Result.fLPowerLawBeta = sqrt(A.fLPowerLawBeta);
	Result.fLPowerLawGamma = sqrt(A.fLPowerLawGamma);
	Result.fBInitMean = sqrt(A.fBInitMean);	
	Result.fBInitSigma = sqrt(A.fBInitSigma);	
	Result.fPInitMean = sqrt(A.fPInitMean);
	Result.fPInitSigma = sqrt(A.fPInitSigma);
	Result.fVKickMod1Mean = sqrt(A.fVKickMod1Mean);
	Result.fVKickMod1Sigma = sqrt(A.fVKickMod1Sigma);
	Result.fVKickMod2Mean = sqrt(A.fVKickMod2Mean);
	Result.fVKickMod2Sigma = sqrt(A.fVKickMod2Sigma);
	Result.fDeathPsi = sqrt(A.fDeathPsi);
	return Result;	
}

SModelParameters operator-(const SModelParameters &A, const SModelParameters &B)
{
	return A + ((-1.)*B);
}

SModelParameters operator+(const SModelParameters &A, const SModelParameters &B)
{
	SModelParameters Result;
	Result.fBDecayDelta    = A.fBDecayDelta    + B.fBDecayDelta;
	Result.fBDecayKappa    = A.fBDecayKappa    + B.fBDecayKappa;
	Result.fLPowerLawAlpha = A.fLPowerLawAlpha + B.fLPowerLawAlpha;
	Result.fLPowerLawBeta  = A.fLPowerLawBeta  + B.fLPowerLawBeta;
	Result.fLPowerLawGamma = A.fLPowerLawGamma + B.fLPowerLawGamma;
	Result.fBInitMean      = A.fBInitMean      + B.fBInitMean;	
	Result.fBInitSigma     = A.fBInitSigma     + B.fBInitSigma;	
	Result.fPInitMean      = A.fPInitMean      + B.fPInitMean;
	Result.fPInitSigma     = A.fPInitSigma     + B.fPInitSigma;
	Result.fVKickMod1Mean  = A.fVKickMod1Mean  + B.fVKickMod1Mean;
	Result.fVKickMod1Sigma = A.fVKickMod1Sigma + B.fVKickMod1Sigma;
	Result.fVKickMod2Mean  = A.fVKickMod2Mean  + B.fVKickMod2Mean;
	Result.fVKickMod2Sigma = A.fVKickMod2Sigma + B.fVKickMod2Sigma;
	Result.fDeathPsi       = A.fDeathPsi       + B.fDeathPsi;

	return Result;	

}

SModelParameters operator*(const SModelParameters &A, const RealType &fB)
{
	return fB*A;
}

SModelParameters operator*(const RealType &fB, const SModelParameters &A)
{
	SModelParameters Result;
	Result.fBDecayDelta = A.fBDecayDelta * fB;
	Result.fBDecayKappa = A.fBDecayKappa * fB;
	Result.fLPowerLawAlpha = A.fLPowerLawAlpha * fB;
	Result.fLPowerLawBeta = A.fLPowerLawBeta * fB;
	Result.fLPowerLawGamma = A.fLPowerLawGamma * fB;
	Result.fBInitMean = A.fBInitMean * fB;	
	Result.fBInitSigma = A.fBInitSigma * fB;	
	Result.fPInitMean = A.fPInitMean * fB;
	Result.fPInitSigma = A.fPInitSigma * fB;
	Result.fVKickMod1Mean = A.fVKickMod1Mean * fB;
	Result.fVKickMod1Sigma = A.fVKickMod1Sigma * fB;
	Result.fVKickMod2Mean = A.fVKickMod2Mean * fB;
	Result.fVKickMod2Sigma = A.fVKickMod2Sigma * fB;
	Result.fDeathPsi = A.fDeathPsi * fB;
	return Result;	
	

}




bool operator<(const SProbabilitySpacePointRadio &A, const SProbabilitySpacePointRadio &B)
{
	if(A.nLogP < B.nLogP) return true;
	if(A.nLogP == B.nLogP)
	{
		if(A.nLogPdot < B.nLogPdot) return true;
		if(A.nLogPdot == B.nLogPdot)
		{
			if(A.nLogS1400 < B.nLogS1400) return true;
		}
	}

	return false;
}


bool operator<(const SProbabilitySpacePointSpace &A, const SProbabilitySpacePointSpace &B)
{
	if(A.nGalL < B.nGalL) return true;
	if(A.nGalL == B.nGalL)
	{
		if(A.nGalKsi < B.nGalKsi) return true;
		if(A.nGalKsi == B.nGalKsi)
		{
			if(A.nDM < B.nDM) return true;
		}
	}
	return false;
}


bool operator<(const SProbabilitySpacePoint4D &A, const SProbabilitySpacePoint4D &B)
{
	if(A.nLogP < B.nLogP) return true;
	if(A.nLogP == B.nLogP)
	{
		if(A.nLogPdot < B.nLogPdot) return true;
		if(A.nLogPdot == B.nLogPdot)
		{
			if(A.nLogS1400 < B.nLogS1400) return true;
			if(A.nLogS1400 == B.nLogS1400)
			{
				if(A.nDM < B.nDM) return true;
				//if(A.nDM == B.nDM)
				//{
				//	if(A.nHEALPix < B.nHEALPix) return true;
				//
				//}
			}
		}
	}
	return false;


}

#include<ConfigContainer.hpp>
SProbabilitySpacePoint4D ComputeSpacePointCoordiantes4D(const SCatalogueEntry PSR)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	SProbabilitySpacePoint4D Result;
 	Result.nLogP      = static_cast<int>(( log10(PSR.fP) - pCfgPtr->Plot.fZeroLogP)*pCfgPtr->Plot.nResLogP/(pCfgPtr->Plot.fMaxLogP-pCfgPtr->Plot.fZeroLogP));
   	Result.nLogPdot   = static_cast<int>(( log10(PSR.fPdot) - pCfgPtr->Plot.fZeroLogPdot)*pCfgPtr->Plot.nResLogPdot/(pCfgPtr->Plot.fMaxLogPdot-pCfgPtr->Plot.fZeroLogPdot));
  	Result.nLogS1400  = static_cast<int>(( log10(PSR.fRadioFlux1400) - pCfgPtr->Plot.fZeroLogS1400)*pCfgPtr->Plot.nResLogS1400/(pCfgPtr->Plot.fMaxLogS1400-pCfgPtr->Plot.fZeroLogS1400));
	Result.nDM     = static_cast<int>( (PSR.fDM * pCfgPtr->Plot.nResDM) / pCfgPtr->Plot.fMaxDM );
	return Result;
}

SProbabilitySpacePoint4D ComputeSpacePointCoordiantes4D(const CPulsar PSR)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	SProbabilitySpacePoint4D Result;
 	Result.nLogP      = static_cast<int>(( log10(PSR.fPeriod) - pCfgPtr->Plot.fZeroLogP)*pCfgPtr->Plot.nResLogP/(pCfgPtr->Plot.fMaxLogP-pCfgPtr->Plot.fZeroLogP));
   	Result.nLogPdot   = static_cast<int>(( log10(PSR.fPeriodDot) - pCfgPtr->Plot.fZeroLogPdot)*pCfgPtr->Plot.nResLogPdot/(pCfgPtr->Plot.fMaxLogPdot-pCfgPtr->Plot.fZeroLogPdot));
  	Result.nLogS1400  = static_cast<int>(( log10(PSR.fS1400) - pCfgPtr->Plot.fZeroLogS1400)*pCfgPtr->Plot.nResLogS1400/(pCfgPtr->Plot.fMaxLogS1400-pCfgPtr->Plot.fZeroLogS1400));
	Result.nDM     = static_cast<int>( (PSR.fDM * pCfgPtr->Plot.nResDM) / pCfgPtr->Plot.fMaxDM );
	return Result;
}

SProbabilitySpacePointRadio ComputeSpacePointCoordiantesRadio(const SCatalogueEntry PSR)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	SProbabilitySpacePointRadio Result;
 	Result.nLogP      = static_cast<int>(( log10(PSR.fP) - pCfgPtr->Plot.fZeroLogP)*pCfgPtr->Plot.nResLogP/(pCfgPtr->Plot.fMaxLogP-pCfgPtr->Plot.fZeroLogP));
   	Result.nLogPdot   = static_cast<int>(( log10(PSR.fPdot) - pCfgPtr->Plot.fZeroLogPdot)*pCfgPtr->Plot.nResLogPdot/(pCfgPtr->Plot.fMaxLogPdot-pCfgPtr->Plot.fZeroLogPdot));
  	Result.nLogS1400  = static_cast<int>(( log10(PSR.fRadioFlux1400) - pCfgPtr->Plot.fZeroLogS1400)*pCfgPtr->Plot.nResLogS1400/(pCfgPtr->Plot.fMaxLogS1400-pCfgPtr->Plot.fZeroLogS1400));


	//cout<<Result.nLogP<<" "<<Result.nLogPdot<<" "<<Result.nLogS1400<<endl;

	return Result;

	
}

SProbabilitySpacePointSpace ComputeSpacePointCoordiantesSpace(const SCatalogueEntry PSR)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	SProbabilitySpacePointSpace Result;
 	Result.nGalL   = static_cast<int>( (PSR.fDegGalL * pCfgPtr->Plot.nResGalL) / 360. );
	Result.nGalKsi = static_cast<int>( (sin((M_PI/180.)*PSR.fDegGalB) +1.) * pCfgPtr->Plot.nResGalKsi /2. );
	Result.nDM     = static_cast<int>( (PSR.fDM * pCfgPtr->Plot.nResDM) / pCfgPtr->Plot.fMaxDM );

//	cout<<Result.nGalL<<" "<<Result.nGalKsi<<" "<<Result.nDM<<endl;
	return Result;
}


SProbabilitySpacePointRadio ComputeSpacePointCoordiantesRadio(const CPulsar PSR)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	SProbabilitySpacePointRadio Result;
 	Result.nLogP      = static_cast<int>(( log10(PSR.fPeriod) - pCfgPtr->Plot.fZeroLogP)*pCfgPtr->Plot.nResLogP/(pCfgPtr->Plot.fMaxLogP-pCfgPtr->Plot.fZeroLogP));
   	Result.nLogPdot   = static_cast<int>(( log10(PSR.fPeriodDot) - pCfgPtr->Plot.fZeroLogPdot)*pCfgPtr->Plot.nResLogPdot/(pCfgPtr->Plot.fMaxLogPdot-pCfgPtr->Plot.fZeroLogPdot));
  	Result.nLogS1400  = static_cast<int>(( log10(PSR.fS1400) - pCfgPtr->Plot.fZeroLogS1400)*pCfgPtr->Plot.nResLogS1400/(pCfgPtr->Plot.fMaxLogS1400-pCfgPtr->Plot.fZeroLogS1400));


	//cout<<Result.nLogP<<" "<<Result.nLogPdot<<" "<<Result.nLogS1400<<endl;

	return Result;



}
 


SProbabilitySpacePointSpace ComputeSpacePointCoordiantesSpace(const CPulsar PSR)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	SProbabilitySpacePointSpace Result;
 	Result.nGalL   = static_cast<int>( (PSR.GalacticPosDeg.fL * pCfgPtr->Plot.nResGalL) / 360. );
	Result.nGalKsi = static_cast<int>( (sin((M_PI/180.)*PSR.GalacticPosDeg.fB) +1.) * pCfgPtr->Plot.nResGalKsi /2. );
	Result.nDM     = static_cast<int>( (PSR.fDM * pCfgPtr->Plot.nResDM) / pCfgPtr->Plot.fMaxDM );

//	cout<<Result.nGalL<<" "<<Result.nGalKsi<<" "<<Result.nDM<<endl;
	return Result;


}

SProbabilitySpacePointRadioSmooth ComputeSpacePointCoordiantesRadioSmooth(const SCatalogueEntry PSR)
{
	SProbabilitySpacePointRadioSmooth Result;
	Result.fLogP = log10(PSR.fP);
	Result.fLogPdot = log10(PSR.fPdot);
	Result.fLogS1400 = log10(PSR.fRadioFlux1400);
	return Result;

}

SProbabilitySpacePointRadioSmooth ComputeSpacePointCoordiantesRadioSmooth(const CPulsar* pPSR)
{
	SProbabilitySpacePointRadioSmooth Result;
	Result.fLogP = log10(pPSR->fPeriod);
	Result.fLogPdot = log10(pPSR->fPeriodDot);
	Result.fLogS1400 = log10(pPSR->fS1400);
	return Result;
}



RealType GaussianishDensityContribution(const SProbabilitySpacePointRadio RefSpacePoint, const SProbabilitySpacePointRadioSmooth SpacePoint, const RealType fSigma)
{
	CConfiguration *pCfgPtr = CConfiguration::GetInstance();
	SProbabilitySpacePointRadioSmooth RefSpacePointSmooth;
 	RefSpacePointSmooth.fLogP = (RefSpacePoint.nLogP*(pCfgPtr->Plot.fMaxLogP-pCfgPtr->Plot.fZeroLogP)/pCfgPtr->Plot.nResLogP) + pCfgPtr->Plot.fZeroLogP;
 	
	RefSpacePointSmooth.fLogPdot = (RefSpacePoint.nLogPdot*(pCfgPtr->Plot.fMaxLogPdot-pCfgPtr->Plot.fZeroLogPdot)/pCfgPtr->Plot.nResLogPdot) + pCfgPtr->Plot.fZeroLogPdot;

 	RefSpacePointSmooth.fLogS1400 = (RefSpacePoint.nLogS1400*(pCfgPtr->Plot.fMaxLogS1400-pCfgPtr->Plot.fZeroLogS1400)/pCfgPtr->Plot.nResLogS1400) + pCfgPtr->Plot.fZeroLogS1400;


    RealType fConst(1./(sqrt(2.*M_PI)*fSigma));	

    RealType fRes = fConst*fConst*fConst* \
            exp(-(RefSpacePointSmooth.fLogP - SpacePoint.fLogP)*(RefSpacePointSmooth.fLogP - SpacePoint.fLogP)/(2.*fSigma*fSigma)) * \
            exp(-(RefSpacePointSmooth.fLogPdot - SpacePoint.fLogPdot)*(RefSpacePointSmooth.fLogPdot - SpacePoint.fLogPdot)/(2.*fSigma*fSigma)) * \
            exp(-(RefSpacePointSmooth.fLogS1400 - SpacePoint.fLogS1400)*(RefSpacePointSmooth.fLogS1400 - SpacePoint.fLogS1400)/(2.*fSigma*fSigma));

//    cout<<fConst<<" "<<fRes<<endl;
//    cout<<"LogP: "<<RefSpacePointSmooth.fLogP<<" "<<SpacePoint.fLogP<<endl;
//    cout<<"LogPdot: "<<RefSpacePointSmooth.fLogPdot<<" "<<SpacePoint.fLogPdot<<endl;
//    cout<<"LogS1400: "<<RefSpacePointSmooth.fLogS1400<<" "<<SpacePoint.fLogS1400<<endl;

    return fRes;
}




