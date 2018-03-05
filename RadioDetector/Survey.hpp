#pragma once
#include<Definitions.hpp>
//FIXME przeniesc do definitions
const RealType gCLight(2.99792458e8);


const int gSurveyNum (10);//bedzie potrzbne w CPulsar dla znalezionych flag
enum ESurveyID {
    Parkes70 = 0,
    ParkesMultiBeam = 1,
    SwinIntermLat = 2,
    SwinExtended = 3,
    Burgay = 4,
    GMRT = 5,
    ParkesMultiBeamALLSKY = 6,
    ParkesMultiBeamPart = 7,
    SKA1Mid = 8,
    SKA1MidLC = 9
};

class CSurvey {
  public:
    // FIXME zmienic nazwy
    static bool SkyCoverageParkes70(const RealType fL, const RealType fB);
    static bool SkyCoverageParkesMB(const RealType fL, const RealType fB);
    static bool SkyCoverageSwinIL(const RealType fL, const RealType fB);
    static bool SkyCoverageSwinExt(const RealType fL, const RealType fB);
    static bool SkyCoverageBurgay(const RealType fL, const RealType fB);
    static bool SkyCoverageGMRT(const RealType fL, const RealType fB);
    static bool SkyCoverageAllSky(const RealType fL, const RealType fB);
    static bool SkyCoveragePartSky(const RealType fL, const RealType fB);
    static bool SkyCoverageSKA1Mid(const RealType fL, const RealType fB);
    static bool SkyCoverageSKA1MidLC(const RealType fL, const RealType fB);
    RealType m_DM0();

    // public:
    ESurveyID ID;                  /// Id-enumerator
    static const int nNameLength = 256;   /// size of the name cstring
    char Name[nNameLength];        /// keep it short
    RealType fTemperatureReceiver; /// Receiver temperature [K]
    RealType fReceiverBandwidth;   /// [MHz] (survey bandwidth)
    int nChannelsNum;              /// Number of channels
    RealType fFrequency;           /// [MHz]
    RealType fSamplingTime;        /// Sampling time [s]
    RealType fGain;                /// Gain of the telescope [K/Jy]
    RealType fIntegrationTime;     /// Integration time [s]
    RealType fDM0;
    /* Pointer to a 'geometry' function. */
    /* should return 1 if (l,b) belongs to this survey, 0 otherwise */
    bool (*IsCovered)(const RealType fL, const RealType fB);

    static bool IsCoveredInSurvey(const RealType fLDeg, const RealType fBDeg,
                                  const ESurveyID SrvID);

    void Init(const ESurveyID SrvID);
    CSurvey(ESurveyID SrvID);
    CSurvey(){};
    ~CSurvey();
	
};



