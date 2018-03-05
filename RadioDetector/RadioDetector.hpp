#pragma once

#include<Definitions.hpp>
#include<Pulsar.hpp>
#include<Survey.hpp>


//FIXME moze aby sie nie chrzanilo zmienic naze na CRadioDetection?
class CRadioDetector {
  protected:
    CSurvey* m_pSurvey;
    int m_nSurveyNum;
    RealType m_SNRTreshold;

    short int *m_pTemperatureSky;
    void m_ReadTemperatureSkyData(const char *Filename);

    RealType m_TauScatter(const CPulsar *pPulsar, const CSurvey* pSurvey) const;
    RealType m_TemperatureSky(const Galactic PosDeg,
                              const CSurvey* pSurvey) const;
    RealType m_TemperatureSkyTab(const RealType fLDeg,
                                 const RealType fBDeg) const;
    RealType m_SNR(const CPulsar *pPulsar, const CSurvey* pSurvey) const;
    RealType m_WEffective(const CPulsar *pPulsar, const CSurvey* pSurvey) const;
    RealType m_WIntrinsic(const CPulsar *pPulsar) const;
    RealType m_SMinimal(const CPulsar *pPulsar,
                        const CSurvey* pSurvey) const; ///<[mJy]
    void m_ComputeDM0(CSurvey* *Survey) const;
    RealType m_Flux(const CPulsar *pPulsar,
                    const CSurvey* pSurvey) const; ///<[mJy] a nie mJy???
    RealType m_Flux400toFreq(const CPulsar *pPulsar,
                             const RealType fFrequency) const; ///<[mJy]
    void CheckVisibility(CPulsar *pPulsar) const;

  public:
    RealType SMinimal(const CPulsar *pPulsar, const ESurveyID ESurvey) const;
    CRadioDetector();
    ~CRadioDetector();
    bool IsVisible(const CPulsar *pPulsar, const ESurveyID ESurvey) const;
    bool CheckParkes(const CPulsar *pPulsar) const;
    void Print4StefansTest(const CPulsar *pPulsar) const;
};

