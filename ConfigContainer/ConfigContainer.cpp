#include<ConfigContainer.hpp>

#include<iostream>
#include<Debug.hpp>
#include<cstring>
#include<sstream>
#include<BasicStar.hpp>
#include<vector>
#include<pugixml.hpp>


using namespace std;


CConfiguration* CConfiguration::GetInstance()
{
	static CConfiguration Instance; // Guaranteed to be destroyed. Instantiated on first use.
	return &Instance;
}


CConfiguration::CConfiguration()
{
	pRan=new CRandomizer;
	Model.nSeed = pRan->m_nSeed;
	BasicConfig();
	Sim.bDynamicsComputed=false;

}

CConfiguration::~CConfiguration()
{
	delete pRan;
	delete pATNF;

	delete pDM;
}

void CConfiguration::Reseed(const unsigned int nSeed, const unsigned nCalls)
{
	pRan->Seed(nSeed, nCalls);
	Model.nSeed = pRan->m_nSeed;

}


//FIXME zaimplementowac czytanie z pliku
void CConfiguration::BasicConfig()
{
    Sim.bStartMCMCFromModel=false;
    Sim.bScaleJumps=true;
    Sim.fJumpScallingFactor=0.75;
    Sim.nScaleEveryJumps=100;

	Sim.fLnLikelihoodSigma = 1e-3;
	Sim.nNumberOfIteration = 5;
	Sim.KickingMethod = HobbsKicking;
	//Sim.ExtremumMethod = KolmogorovSmirnov;
	Sim.ExtremumMethod = Likelihood;
	//Sim.fMaxPulsarAge=1000.;//1e9[Myr]
	Sim.fMaxPulsarAge=50.;//[Myr]
	//Sim.fMaxPulsarAge=25.;//[Myr]
	//Sim.fMaxTimeStep=0.1;//[Myr] Energia całkowita w ruchu w potencjale gal zachowana z dokladnoscia do 1%
	Sim.fMaxTimeStep=0.1;//[Myr] Energia całkowita w ruchu w potencjale gal zachowana z dokladnoscia do 1%
	Sim.nMaxPulsarNum=1000000;
	Sim.nMaxAvailablePulsarNum=1000000;
        Sim.fMaxDist=35.;//[kpc]
        Sim.bObservations=true;
        Sim.bDynamics=true;
        Sim.bPhysics=true;
	Sim.bVerbose=false;
	Sim.bMCMCVerbose=true;
	Sim.bRecycleDynamics=false;
	Sim.bDynamicsComputed=false;
	Sim.bSurveyGeometry=true;
	Sim.bOnlyPowerLaw=true;
	Sim.bStrictLikelihood=false;
	Sim.bATNFDeathLinesCut=false;
	Sim.bCutATNFOutsideArea=true;
	Sim.bComputeCatalogueMaxDist=false;

	Sim.nNumOfChains=1;
	Sim.nNumPerChain=10;


//	Model.fBDecayTimeDeltaMean=25.;//[Myr]
	//Model.fBDecayTimeDeltaMean=4550.59;//[Myr]
//	Model.fBDecayTimeDeltaSigma=0.;//[Myr]


	//TODO Substitude them with best fit 
	Model.fBDecayDelta=10.;
	Model.fBDecayKappa=1.;

	//Set to the same values as Narayan&Ostriker:
	Model.fLPowerLawGamma=5.69e6;
	Model.fLPowerLawAlpha=-1.;
	Model.fLPowerLawBeta=1./3.;

	//rozkłady log-normalne
	Model.fBInitMean=12.65;
	//Model.fBInitSigma=0.55;
	Model.fBInitSigma=1.25;

	Model.fBFinalMean=9.;
	Model.fBFinalSigma=0.5;

	//rozklady liniowe
	Model.fPInitMean=0.300;//[s]
	Model.fPInitSigma=0.150;//[s]

	//DeathArea
	Model.fDeathPsi=0.2;

#if 0
	
	Model.fLPowerLawGamma=1e10;//TODO Co z tym parametrem???
	Model.fLPowerLawAlpha=-1.05;
	Model.fLPowerLawBeta=0.37;

	ModelConstraints.Min.fBDecayDelta=1e-5;
	ModelConstraints.Max.fBDecayDelta=1e5;

	ModelConstraints.Min.fBDecayKappa=0.5;
	ModelConstraints.Max.fBDecayKappa=3.;
	
	ModelConstraints.Min.fLPowerLawAlpha=-1.6;
	ModelConstraints.Max.fLPowerLawAlpha=0.8;

	ModelConstraints.Min.fLPowerLawBeta=0.3;
	ModelConstraints.Max.fLPowerLawBeta=0.7;


#endif
	ModelConstraints.Min.fBDecayDelta=1e-4;//[Myr] -> 100lat
	ModelConstraints.Max.fBDecayDelta=14.4e3;//[Myr] hubbletime
//	ModelConstraints.Sigma.fBDecayDelta=pow(10.,0.05);//[Myr] hubbletime
	

	ModelConstraints.Min.fBDecayKappa=-3.5;
	ModelConstraints.Max.fBDecayKappa=3.5;
//	ModelConstraints.Sigma.fBDecayKappa=0.01;
	
	ModelConstraints.Min.fLPowerLawGamma=1e-3;
	ModelConstraints.Max.fLPowerLawGamma=1e20;
//	ModelConstraints.Sigma.fLPowerLawGamma=pow(10.,0.01);

	ModelConstraints.Min.fLPowerLawAlpha=-2.;
	ModelConstraints.Max.fLPowerLawAlpha=2.;
//	ModelConstraints.Sigma.fLPowerLawAlpha=0.01;

	ModelConstraints.Min.fLPowerLawBeta=-2.;
	ModelConstraints.Max.fLPowerLawBeta=2.;
//	ModelConstraints.Sigma.fLPowerLawBeta=0.01;


	ModelConstraints.Min.fPInitMean=0.150;
	ModelConstraints.Max.fPInitMean=0.450;

	ModelConstraints.Min.fPInitSigma=0.05;
	ModelConstraints.Max.fPInitSigma=0.50;

	ModelConstraints.Min.fBInitMean=10.;
	ModelConstraints.Max.fBInitMean=14.;

	ModelConstraints.Min.fBInitSigma=0.1;
	ModelConstraints.Max.fBInitSigma=5.5;

	ModelConstraints.Min.fBFinalMean=9.;
	ModelConstraints.Max.fBFinalMean=9.;

	ModelConstraints.Min.fBFinalSigma=0.5;
	ModelConstraints.Max.fBFinalSigma=0.5;

	ModelConstraints.Min.fDeathPsi=0.1;
	ModelConstraints.Max.fDeathPsi=5.;


	ModelConstraints.fIntervalFraction=0.1;

	MCMCStepSize.fBDecayDelta=pow(10.,0.1);
	MCMCStepSize.fBDecayKappa=0.1;

	MCMCStepSize.fLPowerLawGamma=pow(10.,0.1);//TODO Co z tym parametrem???
	MCMCStepSize.fLPowerLawAlpha=0.01;
	MCMCStepSize.fLPowerLawBeta=0.01;


	MCMCStepSize.fBInitMean=0.1;
	MCMCStepSize.fBInitSigma=0.1;

	MCMCStepSize.fPInitMean=0.01;
	MCMCStepSize.fPInitSigma=0.01;


	MCMCStepSize.fDeathPsi=0.05;

	//Sim.fIntervalFraction=0.05;
	//Sim.fIntervalScalingFactor=0.95;


	MCMCParameterUse.bLPowerLawGamma = false;
	MCMCParameterUse.bLPowerLawAlpha = false; 
	MCMCParameterUse.bLPowerLawBeta = false;
	MCMCParameterUse.bBDecayDelta = true;
	MCMCParameterUse.bBDecayKappa = false;
	MCMCParameterUse.bBInitMean=true;	
	MCMCParameterUse.bBInitSigma=true;	
	MCMCParameterUse.bPInitMean=true;
	MCMCParameterUse.bPInitSigma=true;
	MCMCParameterUse.bVKickMod1Mean=false;
	MCMCParameterUse.bVKickMod1Sigma=false;
	MCMCParameterUse.bVKickMod2Mean=false;
	MCMCParameterUse.bVKickMod2Sigma=false;
	MCMCParameterUse.bDeathPsi=false;
	MCMCParameterUse.bBFinalMean=false;	
	MCMCParameterUse.bBFinalSigma=false;	

	MCMCParameterConstraint.bLPowerLawGamma = true;
	MCMCParameterConstraint.bLPowerLawAlpha = true; 
	MCMCParameterConstraint.bLPowerLawBeta = true;
	MCMCParameterConstraint.bBDecayDelta = true;
	MCMCParameterConstraint.bBDecayKappa = true;
	MCMCParameterConstraint.bBFinalMean=true;	
	MCMCParameterConstraint.bBFinalSigma=true;	
	MCMCParameterConstraint.bBInitMean=true;	
	MCMCParameterConstraint.bBInitSigma=true;	
	MCMCParameterConstraint.bPInitMean=true;
	MCMCParameterConstraint.bPInitSigma=true;
	MCMCParameterConstraint.bVKickMod1Mean=true;
	MCMCParameterConstraint.bVKickMod1Sigma=true;
	MCMCParameterConstraint.bVKickMod2Mean=true;
	MCMCParameterConstraint.bVKickMod2Sigma=true;
	MCMCParameterConstraint.bDeathPsi=true;

	//0 - Narayan, Ostriker
	//1 - Proszynski, Przybycien
	//2 - PowerLaw
	//3 - Random
	//4 - Erot
	Sim.LuminosityModel=Narayan;

	//0 - Density1x3D 
	//1 - Density2x3D 
	//2 - Density1x3DSparse
	//3 - Density2x3DSparse
	//4 - Density4DSparse
	//5 - Density3x1D
	Sim.ModelDensityContainer=Density1x3D;
	Sim.bGaussianDensityMod=false;
	Sim.fGaussianDensitySigma=1.;

	Sim.bCheck4Relevance = false;

	SurveySelection.nNum=1;///<the number of considered surveys
	SurveySelection.pID= new ESurveyID[SurveySelection.nNum];
	SurveySelection.pID[0]=ParkesMultiBeam;

	Sim.DataFields.Word=0;//zero the datafields' mask

	//set the required datafileds
	Sim.DataFields.Mask.bDegGalL=1;
	Sim.DataFields.Mask.bDegGalB=1;
	Sim.DataFields.Mask.bP=1;
	Sim.DataFields.Mask.bPDot=1;
	Sim.DataFields.Mask.bDM=1;
	Sim.DataFields.Mask.bRadioFlux1400=1;

	Sim.fDensityMinValue=1e-20;

	
	//Plot.nResLogP=128;
	//Plot.nResLogPdot=128;
	//Plot.nResLogS1400=14;
	Plot.nResLogP=30;
	Plot.nResLogPdot=30;
	Plot.nResLogS1400=30;

	Plot.fZeroLogP=-1.5;
	Plot.fZeroLogPdot=-19.;
	Plot.fZeroLogS1400=-2.;
	Plot.fMaxLogP=1.;
	Plot.fMaxLogPdot=-11.;
	Plot.fMaxLogS1400=4.;


	Plot.nResGalL = 32;
	Plot.nResGalKsi = 32;
	Plot.nResDM = 32;

	Plot.fMaxDM = 1500.;



	Sim.bUseCatalogue=true;
	strcpy(Sim.sCatalogueFile, "input/ATNF_MOD1.txt" );
	//pATNF= new CCatalogue();




	strcpy(Sim.sDMDataFile, "input/DMData.bin" );
	
	Sim.bUseDMContainer=false;
	Sim.bUsePrecomputedDMData=false;
	//Sim.bUseDMContainer=true;
	//Sim.bUsePrecomputedDMData=true;
	Sim.nDMDistSlices=350;
	Sim.nHEALPixNSide=128;
	Sim.bComputeDMDiff=false;
	//pDM = new CDMContainer(Sim.sDMDataFile);	
	//pDM = new CDMContainer(Sim.nHEALPixNSide,Sim.nDMDistSlices,Sim.fMaxDist);
	Sim.bVerboseParameterUsage=true;

	Sim.bLnLikelihoodPenalty=false;
	Sim.fLnLikelihoodPenaltyCoefficient=1.;


	Sim.bPhysicalDensity = false;


	Sim.bSameSeedInMCMCModels = false;


	Sim.bLEfficiency = false;
	Sim.bLogLCorr = false;

	Sim.bNoDeathLine = false;

	Sim.fDetectionTreshold = 10.;

}

//TODO dorzucic 

void CConfiguration::PostConfigInit()
{
	if(Sim.bUseCatalogue) pATNF= new CCatalogue();

	if(Sim.bUseDMContainer)
	{
		if(Sim.bUsePrecomputedDMData)
		{
			pDM = new CDMContainer(Sim.sDMDataFile);	
		}
		else
		{
			pDM = new CDMContainer(Sim.nHEALPixNSide,Sim.nDMDistSlices,Sim.fMaxDist);
		}
	}
}


void CConfiguration::Print()
{
	cout<<"Configuration:"<<endl;
	cout<<" MaxPulsarAge: "<<Sim.fMaxPulsarAge<<"[Myr]"<<endl;
	cout<<" MaxTimeStep:  "<<Sim.fMaxTimeStep<<"[Myr]"<<endl;
	cout<<" MaxPulsarNum: "<<Sim.nMaxPulsarNum<<endl;
        cout<<" MaxDistance:  "<<Sim.fMaxDist<<"[kpc]"<<endl;
        cout<<" Dynamics:     "<<Sim.bDynamics<<endl;
        cout<<" Physics:      "<<Sim.bPhysics<<endl;
        cout<<" Observations: "<<Sim.bObservations<<endl;
	cout<<" Seed:         "<<pRan->m_nSeed<<endl;
}
void CConfiguration::ReadXML(const char *FileName)
{
	pugi::xml_document Doc;
	Doc.load_file(FileName);

	cout<<"Loading configuration from: "<<FileName<<endl;
	cout<<"Config name: "<<Doc.child("config").attribute("name").value()<<endl;

	ReadModel(&Doc);

}
void CConfiguration::ReadXMLMCMC(const char *FileName, const bool bVerbose)
{
	pugi::xml_document Doc;
	Doc.load_file(FileName);

	if(bVerbose) cout<<"Loading configuration from: "<<FileName<<endl;
	if(bVerbose) cout<<"Config name: "<<Doc.child("config").attribute("name").value()<<endl;

	ReadModel(&Doc, bVerbose);
	ReadMCMC(&Doc, bVerbose);


}

void CConfiguration::ReadModel(pugi::xml_document *pDoc, const bool bVerbose)
{
	pugi::xml_node Node;


	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bRecycleDynamics");
	Sim.bRecycleDynamics=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bRecycleDynamics "<<Sim.bRecycleDynamics<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bPhysics");
	Sim.bPhysics=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bPhysics "<<Sim.bPhysics<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bObservations");
	Sim.bObservations=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bObservations "<<Sim.bObservations<<endl;
	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bUsePrecomputedDMData");
	Sim.bUsePrecomputedDMData=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bUsePrecomputedDMData "<<Sim.bUsePrecomputedDMData<<endl;


	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bUseDMContainer");
	Sim.bUseDMContainer=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bUseDMContainer "<<Sim.bUseDMContainer<<endl;
	
	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "sDMDataFile");
	strcpy(Sim.sDMDataFile,Node.attribute("value").as_string());
	if(bVerbose) cout<<"sDMDataFile "<<Sim.sDMDataFile<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "sCatalogueFile");
	strcpy(Sim.sCatalogueFile,Node.attribute("value").as_string());
	if(bVerbose) cout<<"sCatalogueFile "<<Sim.sCatalogueFile<<endl;


	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "LuminosityModel");
	Sim.LuminosityModel=static_cast<ELuminosityModel>(Node.attribute("value").as_int());
	if(bVerbose) cout<<"LuminosityModel "<<Sim.LuminosityModel<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fBDecayKappa");
	Model.fBDecayKappa=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fBDecayKappa "<<Model.fBDecayKappa<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fBDecayDelta");
	Model.fBDecayDelta=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fBDecayDelta "<<Model.fBDecayDelta<<endl;
	if(bVerbose) cout<<"Log(Model.fBDecayDelta) "<<log10(Model.fBDecayDelta)<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fLPowerLawAlpha");
	Model.fLPowerLawAlpha=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fLPowerLawAlpha "<<Model.fLPowerLawAlpha<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fLPowerLawBeta");
	Model.fLPowerLawBeta=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fLPowerLawBeta "<<Model.fLPowerLawBeta<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fLPowerLawGamma");
	Model.fLPowerLawGamma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fLPowerLawGamma "<<Model.fLPowerLawGamma<<endl;
	if(bVerbose) cout<<"Model.fLPowerLawGammaLog "<<log10(Model.fLPowerLawGamma)<<endl;


	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "ModelDensityContainer");
	Sim.ModelDensityContainer=static_cast<EModelDensityContainer>(Node.attribute("value").as_int());
	if(bVerbose) cout<<"Sim.DensityContainer "<<Sim.ModelDensityContainer<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bGaussianDensityMod");
	Sim.bGaussianDensityMod=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Sim.bGaussianDensityMod "<<Sim.bGaussianDensityMod<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fGaussianDensitySigma");
	Sim.fGaussianDensitySigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Sim.fGaussianDensitySigma "<<Sim.fGaussianDensitySigma<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fBInitMean");
	Model.fBInitMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fBInitMean "<<Model.fBInitMean<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fBInitSigma");
	Model.fBInitSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fBInitSigma "<<Model.fBInitSigma<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fBFinalMean");
	Model.fBFinalMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fBFinalMean "<<Model.fBFinalMean<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fBFinalSigma");
	Model.fBFinalSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fBFinalSigma "<<Model.fBFinalSigma<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fPInitMean");
	Model.fPInitMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fPInitMean "<<Model.fPInitMean<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fPInitSigma");
	Model.fPInitSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fPInitSigma "<<Model.fPInitSigma<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fDeathPsi");
	Model.fDeathPsi=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Model.fDeathPsi "<<Model.fDeathPsi<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fLnLikelihoodSigma");
	Sim.fLnLikelihoodSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Sim.fLnLikelihoodSigma "<<Sim.fLnLikelihoodSigma<<endl;


	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bLEfficiency");
	Sim.bLEfficiency=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bLEfficiency "<<Sim.bLEfficiency<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bLogLCorr");
	Sim.bLogLCorr=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bLogLCorr "<<Sim.bLogLCorr<<endl;



	Node = pDoc->child("config").child("model").find_child_by_attribute("flag", "varname", "bNoDeathLine");
	Sim.bNoDeathLine=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bNoDeathLine "<<Sim.bNoDeathLine<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "fDetectionTreshold");
	Sim.fDetectionTreshold=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Sim.fDetectionTreshold "<<Sim.fDetectionTreshold<<endl;


	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "nNumberOfIteration");
	Sim.nNumberOfIteration=Node.attribute("value").as_int();
	if(bVerbose) cout<<"nNumberOfIteration "<<Sim.nNumberOfIteration<<endl;

	
	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "nSurveyID");

	//SurveySelection.nNum=1;
	//SurveySelection.pID= new ESurveyID[SurveySelection.nNum];
	//SurveySelection.pID[0]=static_cast<ESurveyID>( Node.attribute("value").as_int() );
	//if(bVerbose) cout<<"nSurveyID "<<SurveySelection.pID[0]<<endl;

	Node = pDoc->child("config").child("model").find_child_by_attribute("parameter", "varname", "sSurveysID");
	std::stringstream sSurveysID(Node.attribute("value").as_string());

	int nTempVal(-1);	
	std::vector<ESurveyID> TempVec;
	
	while(sSurveysID>>nTempVal) 
	{
		TempVec.push_back(static_cast<ESurveyID>(nTempVal));
	}

	if(SurveySelection.pID) 
	{
		delete [] SurveySelection.pID;
	}
	SurveySelection.pID= new ESurveyID[TempVec.size()];
	
	for(int i(0);i<TempVec.size();i++)
	{
		SurveySelection.pID[i]=TempVec[i];
		if(bVerbose) cout<<"SurveyID: "<<i<<" "<<SurveySelection.pID[i]<<endl;
	}	

	SurveySelection.nNum=TempVec.size();

}

void CConfiguration::ReadMCMC(pugi::xml_document *pDoc, const bool bVerbose)
{

	
	pugi::xml_node Node;

	if(bVerbose) cout<<"============================="<<endl;
	if(bVerbose) cout<<"Chain geometry per simulation"<<endl;	
	if(bVerbose) cout<<"============================="<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "nNumOfChains");
	Sim.nNumOfChains=Node.attribute("value").as_int();
	if(bVerbose) cout<<"nNumOfChains "<<Sim.nNumOfChains<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "nNumPerChain");
	Sim.nNumPerChain=Node.attribute("value").as_int();
	if(bVerbose) cout<<"nNumPerChain "<<Sim.nNumPerChain<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "bLnLikelihoodPenalty");
	Sim.bLnLikelihoodPenalty=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bLnLikelihoodPenalty "<<Sim.bLnLikelihoodPenalty<<endl;

	
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "bScaleJumps");
	Sim.bScaleJumps=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bScaleJumps "<<Sim.bScaleJumps<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "nScaleEveryJumps");
	Sim.nScaleEveryJumps=Node.attribute("value").as_int();
	if(bVerbose) cout<<"nScaleEveryJumps "<<Sim.nScaleEveryJumps<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "fJumpScallingFactor");
	Sim.fJumpScallingFactor=Node.attribute("value").as_double();
	if(bVerbose) cout<<"fJumpScallingFactor "<<Sim.fJumpScallingFactor<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "bStartMCMCFromModel");
	Sim.bStartMCMCFromModel=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"bStartMCMCFromModel "<<Sim.bStartMCMCFromModel<<endl;









	if(bVerbose) cout<<"=================="<<endl;
	if(bVerbose) cout<<"Model flags' usage"<<endl;
	if(bVerbose) cout<<"=================="<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bLPowerLawGamma");
	MCMCParameterUse.bLPowerLawGamma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bLPowerLawGamma "<<MCMCParameterUse.bLPowerLawGamma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bLPowerLawAlpha");
	MCMCParameterUse.bLPowerLawAlpha=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bLPowerLawAlpha "<<MCMCParameterUse.bLPowerLawAlpha<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bLPowerLawBeta");
	MCMCParameterUse.bLPowerLawBeta=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bLPowerLawBeta "<<MCMCParameterUse.bLPowerLawBeta<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bBDecayDelta");
	MCMCParameterUse.bBDecayDelta=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bBDecayDelta "<<MCMCParameterUse.bBDecayDelta<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bBDecayKappa");
	MCMCParameterUse.bBDecayKappa=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bBDecayKappa "<<MCMCParameterUse.bBDecayKappa<<endl;



	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bBInitMean");
	MCMCParameterUse.bBInitMean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bBInitMean "<<MCMCParameterUse.bBInitMean<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bBInitSigma");
	MCMCParameterUse.bBInitSigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bBInitSigma "<<MCMCParameterUse.bBInitSigma<<endl;


	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bBFinalMean");
	MCMCParameterUse.bBFinalMean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bBFinalMean "<<MCMCParameterUse.bBFinalMean<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bBFinalSigma");
	MCMCParameterUse.bBFinalSigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bBInitSigma "<<MCMCParameterUse.bBInitSigma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bPInitMean");
	MCMCParameterUse.bPInitMean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bPInitMean "<<MCMCParameterUse.bPInitMean<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bPInitSigma");
	MCMCParameterUse.bPInitSigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bPInitSigma "<<MCMCParameterUse.bPInitSigma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bVKickMod1Mean");
	MCMCParameterUse.bVKickMod1Mean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bVKickMod1Mean "<<MCMCParameterUse.bVKickMod1Mean<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bVKickMod1Sigma");
	MCMCParameterUse.bVKickMod1Sigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bVKickMod1Sigma "<<MCMCParameterUse.bVKickMod1Sigma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bVKickMod2Mean");
	MCMCParameterUse.bVKickMod2Mean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bVKickMod2Mean "<<MCMCParameterUse.bVKickMod2Mean<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bVKickMod2Sigma");
	MCMCParameterUse.bVKickMod2Sigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bVKickMod2Sigma "<<MCMCParameterUse.bVKickMod2Sigma<<endl;


	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Use.bDeathPsi");
	MCMCParameterUse.bDeathPsi=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Use.bDeathPsi "<<MCMCParameterUse.bDeathPsi<<endl;










	if(bVerbose) cout<<"============================"<<endl;
	if(bVerbose) cout<<"Model Parameters' constraint"<<endl;
	if(bVerbose) cout<<"============================"<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bLPowerLawGamma");
	MCMCParameterConstraint.bLPowerLawGamma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bLPowerLawGamma "<<MCMCParameterConstraint.bLPowerLawGamma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bLPowerLawAlpha");
	MCMCParameterConstraint.bLPowerLawAlpha=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bLPowerLawAlpha "<<MCMCParameterConstraint.bLPowerLawAlpha<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bLPowerLawBeta");
	MCMCParameterConstraint.bLPowerLawBeta=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bLPowerLawBeta "<<MCMCParameterConstraint.bLPowerLawBeta<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bBDecayDelta");
	MCMCParameterConstraint.bBDecayDelta=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bBDecayDelta "<<MCMCParameterConstraint.bBDecayDelta<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bBDecayKappa");
	MCMCParameterConstraint.bBDecayKappa=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bBDecayKappa "<<MCMCParameterConstraint.bBDecayKappa<<endl;



	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bBInitMean");
	MCMCParameterConstraint.bBInitMean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bBInitMean "<<MCMCParameterConstraint.bBInitMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bBInitSigma");
	MCMCParameterConstraint.bBInitSigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bBInitSigma "<<MCMCParameterConstraint.bBInitSigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bBFinalMean");
	MCMCParameterConstraint.bBFinalMean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bBFinalMean "<<MCMCParameterConstraint.bBFinalMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bBFinalSigma");
	MCMCParameterConstraint.bBFinalSigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bBFinalSigma "<<MCMCParameterConstraint.bBFinalSigma<<endl;



	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bPInitMean");
	MCMCParameterConstraint.bPInitMean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bPInitMean "<<MCMCParameterConstraint.bPInitMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bPInitSigma");
	MCMCParameterConstraint.bPInitSigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bPInitSigma "<<MCMCParameterConstraint.bPInitSigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bVKickMod1Mean");
	MCMCParameterConstraint.bVKickMod1Mean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bVKickMod1Mean "<<MCMCParameterConstraint.bVKickMod1Mean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bVKickMod1Sigma");
	MCMCParameterConstraint.bVKickMod1Sigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bVKickMod1Sigma "<<MCMCParameterConstraint.bVKickMod1Sigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bVKickMod2Mean");
	MCMCParameterConstraint.bVKickMod2Mean=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bVKickMod2Mean "<<MCMCParameterConstraint.bVKickMod2Mean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bVKickMod2Sigma");
	MCMCParameterConstraint.bVKickMod2Sigma=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bVKickMod2Sigma "<<MCMCParameterConstraint.bVKickMod2Sigma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("flag", "varname", "Constraint.bDeathPsi");
	MCMCParameterConstraint.bDeathPsi=Node.attribute("value").as_bool();
	if(bVerbose) cout<<"Constraint.bDeathPsi "<<MCMCParameterConstraint.bDeathPsi<<endl;








	if(bVerbose) cout<<"==================="<<endl;
	if(bVerbose) cout<<"Constraints values:"<<endl;
	if(bVerbose) cout<<"==================="<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fBDecayDelta");
	ModelConstraints.Min.fBDecayDelta=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fBDecayDelta "<<ModelConstraints.Min.fBDecayDelta<<endl;
	if(bVerbose) cout<<"Min.fBDecayDeltaLog "<<log10(ModelConstraints.Min.fBDecayDelta)<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fBDecayDelta");
	ModelConstraints.Max.fBDecayDelta=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fBDecayDelta "<<ModelConstraints.Max.fBDecayDelta<<endl;
	if(bVerbose) cout<<"Max.fBDecayDeltaLog "<<log10(ModelConstraints.Max.fBDecayDelta)<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fBDecayKappa");
	ModelConstraints.Min.fBDecayKappa=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fBDecayKappa "<<ModelConstraints.Min.fBDecayKappa<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fBDecayKappa");
	ModelConstraints.Max.fBDecayKappa=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fBDecayKappa "<<ModelConstraints.Max.fBDecayKappa<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fLPowerLawGamma");
	ModelConstraints.Min.fLPowerLawGamma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fLPowerLawGamma "<<ModelConstraints.Min.fLPowerLawGamma<<endl;
	if(bVerbose) cout<<"Min.fLPowerLawGammaLog "<<log10(ModelConstraints.Min.fLPowerLawGamma)<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fLPowerLawGamma");
	ModelConstraints.Max.fLPowerLawGamma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fLPowerLawGamma "<<ModelConstraints.Max.fLPowerLawGamma<<endl;
	if(bVerbose) cout<<"Max.fLPowerLawGammaLog "<<log10(ModelConstraints.Max.fLPowerLawGamma)<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fLPowerLawAlpha");
	ModelConstraints.Min.fLPowerLawAlpha=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fLPowerLawAlpha "<<ModelConstraints.Min.fLPowerLawAlpha<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fLPowerLawAlpha");
	ModelConstraints.Max.fLPowerLawAlpha=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fLPowerLawAlpha "<<ModelConstraints.Max.fLPowerLawAlpha<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fLPowerLawBeta");
	ModelConstraints.Min.fLPowerLawBeta=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fLPowerLawBeta "<<ModelConstraints.Min.fLPowerLawBeta<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fLPowerLawBeta");
	ModelConstraints.Max.fLPowerLawBeta=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fLPowerLawBeta "<<ModelConstraints.Max.fLPowerLawBeta<<endl;



	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fBInitMean");
	ModelConstraints.Max.fBInitMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fBInitMean "<<ModelConstraints.Max.fBInitMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fBInitSigma");
	ModelConstraints.Max.fBInitSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fBInitSigma "<<ModelConstraints.Max.fBInitSigma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fBFinalMean");
	ModelConstraints.Max.fBFinalMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fBFinalMean "<<ModelConstraints.Max.fBFinalMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fBFinalSigma");
	ModelConstraints.Max.fBFinalSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fBFinalSigma "<<ModelConstraints.Max.fBFinalSigma<<endl;


	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fPInitMean");
	ModelConstraints.Max.fPInitMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fPInitMean "<<ModelConstraints.Max.fPInitMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fPInitSigma");
	ModelConstraints.Max.fPInitSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fPInitSigma "<<ModelConstraints.Max.fPInitSigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fVKickMod1Mean");
	ModelConstraints.Max.fVKickMod1Mean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fVKickMod1Mean "<<ModelConstraints.Max.fVKickMod1Mean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fVKickMod1Sigma");
	ModelConstraints.Max.fVKickMod1Sigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fVKickMod1Sigma "<<ModelConstraints.Max.fVKickMod1Sigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fVKickMod2Mean");
	ModelConstraints.Max.fVKickMod2Mean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fVKickMod2Mean "<<ModelConstraints.Max.fVKickMod2Mean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fVKickMod2Sigma");



	ModelConstraints.Min.fVKickMod2Sigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fVKickMod2Sigma "<<ModelConstraints.Min.fVKickMod2Sigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fBInitMean");
	ModelConstraints.Min.fBInitMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fBInitMean "<<ModelConstraints.Min.fBInitMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fBInitSigma");
	ModelConstraints.Min.fBInitSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fBInitSigma "<<ModelConstraints.Min.fBInitSigma<<endl;


	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fBFinalMean");
	ModelConstraints.Min.fBFinalMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fBFinalMean "<<ModelConstraints.Min.fBFinalMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fBFinalSigma");
	ModelConstraints.Min.fBFinalSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fBFinalSigma "<<ModelConstraints.Min.fBFinalSigma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fPInitMean");
	ModelConstraints.Min.fPInitMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fPInitMean "<<ModelConstraints.Min.fPInitMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fPInitSigma");
	ModelConstraints.Min.fPInitSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fPInitSigma "<<ModelConstraints.Min.fPInitSigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fVKickMod1Mean");
	ModelConstraints.Min.fVKickMod1Mean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fVKickMod1Mean "<<ModelConstraints.Min.fVKickMod1Mean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fVKickMod1Sigma");
	ModelConstraints.Min.fVKickMod1Sigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fVKickMod1Sigma "<<ModelConstraints.Min.fVKickMod1Sigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fVKickMod2Mean");
	ModelConstraints.Min.fVKickMod2Mean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fVKickMod2Mean "<<ModelConstraints.Min.fVKickMod2Mean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fVKickMod2Sigma");
	ModelConstraints.Min.fVKickMod2Sigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fVKickMod2Sigma "<<ModelConstraints.Min.fVKickMod2Sigma<<endl;


	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Min.fDeathPsi");
	ModelConstraints.Min.fDeathPsi=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Min.fDeathPsi "<<ModelConstraints.Min.fDeathPsi<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "Max.fDeathPsi");
	ModelConstraints.Max.fDeathPsi=Node.attribute("value").as_double();
	if(bVerbose) cout<<"Max.fDeathPsi "<<ModelConstraints.Max.fDeathPsi<<endl;




	if(bVerbose) cout<<"========="<<endl;
	if(bVerbose) cout<<"Step size"<<endl;
	if(bVerbose) cout<<"========="<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fBDecayDelta");
	MCMCStepSize.fBDecayDelta=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fBDecayDelta "<<MCMCStepSize.fBDecayDelta<<endl;
	if(bVerbose) cout<<"MCMCStepSize.fBDecayDeltaLog "<<log10(MCMCStepSize.fBDecayDelta)<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fBDecayKappa");
	MCMCStepSize.fBDecayKappa=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fBDecayKappa "<<MCMCStepSize.fBDecayKappa<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fLPowerLawGamma");
	MCMCStepSize.fLPowerLawGamma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fLPowerLawGamma "<<MCMCStepSize.fLPowerLawGamma<<endl;
	if(bVerbose) cout<<"MCMCStepSize.fLPowerLawGammaLog "<<log10(MCMCStepSize.fLPowerLawGamma)<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fLPowerLawAlpha");
	MCMCStepSize.fLPowerLawAlpha=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fLPowerLawAlpha "<<MCMCStepSize.fLPowerLawAlpha<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fLPowerLawBeta");
	MCMCStepSize.fLPowerLawBeta=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fLPowerLawBeta "<<MCMCStepSize.fLPowerLawBeta<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fBInitMean");
	MCMCStepSize.fBInitMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fBInitMean "<<MCMCStepSize.fBInitMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fBInitSigma");
	MCMCStepSize.fBInitSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fBInitSigma "<<MCMCStepSize.fBInitSigma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fBFinalMean");
	MCMCStepSize.fBFinalMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fBFinalMean "<<MCMCStepSize.fBFinalMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fBFinalSigma");
	MCMCStepSize.fBFinalSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fBFinalSigma "<<MCMCStepSize.fBFinalSigma<<endl;

	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fPInitMean");
	MCMCStepSize.fPInitMean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fPInitMean "<<MCMCStepSize.fPInitMean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fPInitSigma");
	MCMCStepSize.fPInitSigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fPInitSigma "<<MCMCStepSize.fPInitSigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fVKickMod1Mean");
	MCMCStepSize.fVKickMod1Mean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fVKickMod1Mean "<<MCMCStepSize.fVKickMod1Mean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fVKickMod1Sigma");
	MCMCStepSize.fVKickMod1Sigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fVKickMod1Sigma "<<MCMCStepSize.fVKickMod1Sigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fVKickMod2Mean");
	MCMCStepSize.fVKickMod2Mean=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fVKickMod2Mean "<<MCMCStepSize.fVKickMod2Mean<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fVKickMod2Sigma");
	MCMCStepSize.fVKickMod2Sigma=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fVKickMod2Sigma "<<MCMCStepSize.fVKickMod2Sigma<<endl;
	Node = pDoc->child("config").child("MCMC").find_child_by_attribute("parameter", "varname", "MCMCStepSize.fDeathPsi");
	MCMCStepSize.fDeathPsi=Node.attribute("value").as_double();
	if(bVerbose) cout<<"MCMCStepSize.fDeathPsi "<<MCMCStepSize.fDeathPsi<<endl;

}
//TODO dopisac zapis wielkosci Log
void CConfiguration::WriteXML(const char *FileName)
{
	pugi::xml_document XMLDoc;
	pugi::xml_node Config = XMLDoc.append_child("config");
	Config.append_attribute("name") = "Basic Configuration";


	XMLDoc.save_file(FileName);

}

void CConfiguration::AppendModel(pugi::xml_node *pConfig)
{
	pugi::xml_node Simulation = pConfig->append_child("model");
	Simulation.append_attribute("comment") = "Default/best model parameters";

	pugi::xml_node Node;

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "ModelDensityContainer";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "int";
	Node.append_attribute("value") =  Sim.ModelDensityContainer;
	Node.append_attribute("comment") = "Density1x3D=0, Density2x3D=1, Density1x3DSparse=2, Density2x3DSparse=3, Density4DSparse=4, Density3x1D=5";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bGaussianDensityMod";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = Sim.bGaussianDensityMod ;
	Node.append_attribute("comment") = "If the Density1x3D should use the Gaussian patch.";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fGaussianDensitySigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = Sim.fGaussianDensitySigma;
	Node.append_attribute("comment") = "Sigma for Gaussian patch.";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fLPowerLawGamma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = Model.fLPowerLawGamma;
	Node.append_attribute("comment") = "PowerLaw scale parameter [depends on alpha and beta]";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fLPowerLawAlpha";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fLPowerLawAlpha;
	Node.append_attribute("comment") = "PowerLaw spin period exponent";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fLPowerLawBeta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = Model.fLPowerLawBeta;
	Node.append_attribute("comment") = "PowerLaw spin period derivative*1e-15 exponent";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fBDecayDelta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = Model.fBDecayDelta;
	Node.append_attribute("comment") = "Magnetic field decay scale [Myr]";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fBDecayKappa";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fBDecayKappa;
	Node.append_attribute("comment") = "Magnetic field decay exponent";


	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fBInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fBInitMean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fBInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fBInitSigma;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fBFinalMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fBFinalMean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fBFinalSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fBFinalSigma;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fPInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fPInitMean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fPInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fPInitSigma;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fVKickMod1Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fVKickMod1Mean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fVKickMod1Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fVKickMod1Sigma;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fVKickMod2Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fVKickMod2Mean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fVKickMod2Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fVKickMod2Sigma;
	Node.append_attribute("comment") = "";


	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fDeathPsi";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Model.fDeathPsi;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bRecycleDynamics";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = Sim.bRecycleDynamics ;
	Node.append_attribute("comment") = "If model should use geometrical population or compute it each time";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bPhysics";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = Sim.bPhysics ;
	Node.append_attribute("comment") = "If physical model should be compted";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bObservations";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = Sim.bObservations ;
	Node.append_attribute("comment") = "If model should be compared to a survey";


	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "LuminosityModel";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "int";
	Node.append_attribute("value") =  Sim.LuminosityModel;
	Node.append_attribute("comment") = "Narayan=0, Proszynski=1, PowerLaw=2, Random=3, Erot=4";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bUsePrecomputedDMData";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") =  Sim.bUsePrecomputedDMData;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bUseDMContainer";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") =  Sim.bUseDMContainer;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "sDMDataFile";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "string";
	Node.append_attribute("value") =  Sim.sDMDataFile;
	Node.append_attribute("comment") = "";


	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "sCatalogueFile";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "string";
	Node.append_attribute("value") =  Sim.sCatalogueFile;
	Node.append_attribute("comment") = "";



	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bLEfficiency";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") =  Sim.bLEfficiency;
	Node.append_attribute("comment") = "";


	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bLogLCorr";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") =  Sim.bLogLCorr;
	Node.append_attribute("comment") = "";


	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bNoDeathLine";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") =  Sim.bNoDeathLine;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fDetectionTreshold";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Sim.fDetectionTreshold;
	Node.append_attribute("comment") = "The orig. value is 10.";



	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fLnLikelihoodSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") =  Sim.fLnLikelihoodSigma;
	Node.append_attribute("comment") = "Default 1e-3";





	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "nNumberOfIteration";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "int";
	Node.append_attribute("value") =  Sim.nNumberOfIteration;
	Node.append_attribute("comment") = "The number of iterations. One position is represented for n-PSRs. The default value is 5.";







	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "sSurveysID";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "string";

	std::stringstream sSurveysID;

	for(int i(0);i<SurveySelection.nNum;i++)
	{
		sSurveysID<<SurveySelection.pID[i]<<" ";
	}

	Node.append_attribute("value") =  sSurveysID.str().c_str();
	Node.append_attribute("comment") = "The Survey ID. Parkes70=0, ParkesMultiBeam=1, SwinIntermLat=2, SwinExtended=3, Burgay=4, GMRT=5, ParkesMultiBeamALLSKY=6, ParkesMultiBeamPart=7, SKA1Mid=8, SKA1MidLC=9. The default value is 1 (ParkesMB).";

}

void CConfiguration::AppendMCMC(pugi::xml_node *pConfig)
{
	pugi::xml_node Simulation = pConfig->append_child("MCMC");
	Simulation.append_attribute("comment") = "Markov chain Monte Carlo parametrs";

	pugi::xml_node Node;

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "nNumOfChains";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "int";
	Node.append_attribute("value") = Sim.nNumOfChains;
	Node.append_attribute("comment") = "Number of chains to compute";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "nNumPerChain";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "int";
	Node.append_attribute("value") = Sim.nNumPerChain;
	Node.append_attribute("comment") = "Number of links per chain";

    Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "nScaleEveryJumps";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "int";
	Node.append_attribute("value") = Sim.nScaleEveryJumps;
	Node.append_attribute("comment") = "The scaling factor to be applied every n-th jump";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bScaleJumps";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") =  Sim.bScaleJumps;
	Node.append_attribute("comment") = "Scale jumps sigmas with the scaling factor.";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "fJumpScallingFactor";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = Sim.fJumpScallingFactor;
	Node.append_attribute("comment") = "The scaling factor for jumps.";

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bStartMCMCFromModel";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") =  Sim.bStartMCMCFromModel;
	Node.append_attribute("comment") = "Don't use the uniform initial distr. start from the model.";

#if 0
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = " ";
	Node.append_attribute("value") =  ;
	Node.append_attribute("comment") = " ";
#endif

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "bLnLikelihoodPenalty";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") =  Sim.bLnLikelihoodPenalty;
	Node.append_attribute("comment") = "Use the penalization druign LnLikelihood computation";





	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fBDecayDelta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fBDecayDelta;
	Node.append_attribute("comment") = "Min of decay scale [Myr]";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fBDecayDelta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fBDecayDelta;
	Node.append_attribute("comment") = "Max of decay scale [Myr]";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fBDecayKappa";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fBDecayKappa;
	Node.append_attribute("comment") = "Min of decay exponent";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fBDecayKappa";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fBDecayKappa;
	Node.append_attribute("comment") = "Max of decay exponent";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fLPowerLawGamma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fLPowerLawGamma;
	Node.append_attribute("comment") = "Min of PowerLaw scale parameter [depends on alpha and beta]";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fLPowerLawGamma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fLPowerLawGamma;
	Node.append_attribute("comment") = "Max of PowerLaw scale parameter [depends on alpha and beta]";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fLPowerLawAlpha";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fLPowerLawAlpha;
	Node.append_attribute("comment") = "Min of PowerLaw spin period exponent";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fLPowerLawAlpha";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fLPowerLawAlpha;
	Node.append_attribute("comment") = "Max of PowerLaw spin period exponent";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fLPowerLawBeta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fLPowerLawBeta;
	Node.append_attribute("comment") = "Min of PowerLaw spin period derivative*1e-15 exponent";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fLPowerLawBeta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fLPowerLawBeta;
	Node.append_attribute("comment") = "Max of PowerLaw spin period derivative*1e-15 exponent";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fBInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fBInitMean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fBInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fBInitMean;
	Node.append_attribute("comment") = "";
	
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fBFinalMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fBFinalMean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fBFinalMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fBFinalMean;
	Node.append_attribute("comment") = "";


	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fBInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fBInitSigma;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fBInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fBInitSigma;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fPInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fPInitMean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fPInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fPInitMean;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fPInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fPInitSigma;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fPInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fPInitSigma;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fVKickMod1Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fVKickMod1Mean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fVKickMod1Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fVKickMod1Mean;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fVKickMod1Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fVKickMod1Sigma;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fVKickMod1Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fVKickMod1Sigma;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fVKickMod2Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fVKickMod2Mean;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fVKickMod2Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fVKickMod2Mean;
	Node.append_attribute("comment") = "";

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fVKickMod2Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fVKickMod2Sigma;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fVKickMod2Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fVKickMod2Sigma;
	Node.append_attribute("comment") = "";


	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Min.fDeathPsi";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Min.fDeathPsi;
	Node.append_attribute("comment") = "";
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "Max.fDeathPsi";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = ModelConstraints.Max.fDeathPsi;
	Node.append_attribute("comment") = "";




#if 0
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = " ";
	Node.append_attribute("value") =  ;
	Node.append_attribute("comment") = " ";
#endif

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bLPowerLawGamma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bLPowerLawGamma ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bLPowerLawAlpha";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bLPowerLawAlpha ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bLPowerLawBeta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bLPowerLawBeta ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bBDecayDelta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bBDecayDelta ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bBDecayKappa";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bBDecayKappa ;



	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bBInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bBInitMean ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bBInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bBInitSigma ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bBFinalMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bBFinalMean ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bBFinalSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bBFinalSigma ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bPInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bPInitMean ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bPInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bPInitSigma ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bVKickMod1Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bVKickMod1Mean ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bVKickMod1Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bVKickMod1Sigma ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bVKickMod2Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bVKickMod2Mean ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bVKickMod2Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bVKickMod2Sigma ;


	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Use.bDeathPsi";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterUse.bDeathPsi ;











	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bLPowerLawGamma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bLPowerLawGamma ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bLPowerLawAlpha";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bLPowerLawAlpha ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bLPowerLawBeta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bLPowerLawBeta ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bBDecayDelta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bBDecayDelta ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bBDecayKappa";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bBDecayKappa ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bBInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bBInitMean ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bBInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bBInitSigma ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bBFinalMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bBFinalMean ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bBFinalSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bBFinalSigma ;

	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bPInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bPInitMean ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bPInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bPInitSigma ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bVKickMod1Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bVKickMod1Mean ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bVKickMod1Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bVKickMod1Sigma ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bVKickMod2Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bVKickMod2Mean ;
	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bVKickMod2Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bVKickMod2Sigma ;


	Node = Simulation.append_child("flag");
	Node.append_attribute("varname") = "Constraint.bDeathPsi";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "bool";
	Node.append_attribute("value") = MCMCParameterConstraint.bDeathPsi ;




	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fBDecayDelta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fBDecayDelta;
	
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fBDecayKappa";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fBDecayKappa;



	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fLPowerLawGamma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fLPowerLawGamma;

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fLPowerLawAlpha";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fLPowerLawAlpha;

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fLPowerLawBeta";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fLPowerLawBeta;

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fBInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fBInitMean;
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fBInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fBInitSigma;

	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fBFinalMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fBFinalMean;
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fBFinalSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fBFinalSigma;
	Node = Simulation.append_child("parameter");

	Node.append_attribute("varname") = "MCMCStepSize.fPInitMean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fPInitMean;
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fPInitSigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fPInitSigma;
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fVKickMod1Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fVKickMod1Mean;
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fVKickMod1Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fVKickMod1Sigma;
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fVKickMod2Mean";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fVKickMod2Mean;
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fVKickMod2Sigma";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fVKickMod2Sigma;
	Node = Simulation.append_child("parameter");
	Node.append_attribute("varname") = "MCMCStepSize.fDeathPsi";
	Node.insert_attribute_after("type", Node.attribute("varname")) = "double";
	Node.append_attribute("value") = MCMCStepSize.fDeathPsi;


}


void CConfiguration::WriteXMLMCMC(const char *FileName)
{
	pugi::xml_document XMLDoc;
	pugi::xml_node Config = XMLDoc.append_child("config");
	Config.append_attribute("name") = "Single node configuration";


	AppendMCMC(&Config);

	AppendModel(&Config);

	XMLDoc.save_file(FileName);

}
