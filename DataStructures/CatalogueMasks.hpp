#pragma once
#include<cstdint>
#include<string>



//FIXME do zintegrowania z obserwerem
//FIXME przeniesc do datatypes
union SSurveys
{
	struct SMask
	{	
		bool bParkes70             :1;
		bool bParkesMultiBeam      :1;
		bool bSwinIntermLat        :1;
		bool bSwinExtended         :1;
		bool bBurgay               :1;
		bool bGMRT                 :1;
		bool bParkesMultiBeamAllSky:1;
		bool bParkesMultiBeamPart  :1;
	} Mask;
	unsigned char Word;
};


union SDataFields
{
	struct SMask
	{
		bool bDegGalL            :1;
		bool bDegGalB            :1;
		bool bP                  :1;
		bool bPDot               :1;
		bool bDM                 :1;
		bool bRM                 :1;
		bool bDist               :1;
		bool bDistDM             :1;
		bool bRadioFlux400       :1;
		bool bRadioFlux1400      :1;
		bool bAge                :1;
		bool bBSurface           :1;
		unsigned char Fill       :4;
	} Mask;
	uint16_t Word;
};


//to musiało gdzieś wylądować...
struct SCatalogueEntry
{
	SSurveys Surveys;
	SDataFields DataFields;
	std::string Name;
	int nSurveys;//do zastapienia SuveryMask
	RealType fDegGalL;
	RealType fDegGalB;
	RealType fP;
	RealType fPdot;
	RealType fDM;
	RealType fRM;
	RealType fDist;
	RealType fDistDM;
	RealType fRadioFlux400;//[mJy]
	RealType fRadioFlux1400;//[mJy]
	RealType fAge;
	RealType fBSurface;
};


