#pragma once

#include<ConfigContainer.hpp>
#include<CatalogueMasks.hpp>
#include<ModelParameters.hpp>
#include<Definitions.hpp>
#include<map>
#include<Statistics.hpp>
using namespace std;


class CCatalogue
{
private:
	CConfiguration *m_pCfgPtr;
	
public:


	void ReadTxt(const char * FileName);
	void ReadBin(const char * FileName);
	void WriteBin(const char * FileName);	

	SSurveys Check4Surveys(const string Input) const;

	int nTotalNum;
	SCatalogueEntry *pPulsar;
    bool *pbUsed;

	CCatalogue(const char * FileName);
	CCatalogue(); ///Assumes defult catalogue file
	~CCatalogue();

	void Print();


};





SSurveys SurveyID2Word(ESurveyID ID);

class CObservations
{
private:
	void m_PrintUsedEntries(const bool *pMaska);
	int nUsedATNFPulsars;
	CConfiguration *m_pCfgPtr;
	void m_AddToMap(const SProbabilitySpacePointRadio SpacePoint);
	void m_AddToMap(const SProbabilitySpacePointSpace SpacePoint);
	void m_AddToMap(const SProbabilitySpacePoint4D Point);
	bool m_AboveDeathLines(const SCatalogueEntry Pulsar) const; //returns true if above death line
	bool m_InsideArea(const SCatalogueEntry Pulsar, const ESurveyID ID) const;

	void m_AddToCM(const SCatalogueEntry Pulsar);
	void m_ComputeCM(const int nTotalMass);

	void NormContainers();
public:
	CDensityVector Vector;
	RealType fMaxDist;
	SSurveys	SurveyMask;
	ESurveyID	SurveyID;	
	//SCatalogueEntry **ppPSRObs;
	CSurvey Survey;
	std::map<SProbabilitySpacePointRadio, RealType> MapRadio;
	std::map<SProbabilitySpacePointSpace, RealType> MapSpace;

	std::map<SProbabilitySpacePoint4D, RealType> Map4D;

	RealType ComputeMaxDist(const ESurveyID ID = ParkesMultiBeam);
	void Write(const char* Prefix = NULL);
	CObservations(const ESurveyID SurveyID);
	~CObservations(){};
	
	SCenterOfMassRadio CMRadio;
	SCenterOfMass4D    CM4D;

	RealType DistanceToCM(const CPulsar PSR);
};






class CTestObservations: public CObservations{
public:
	CTestObservations():CObservations(ParkesMultiBeam){};
	~CTestObservations(){};

	void Test();

};
