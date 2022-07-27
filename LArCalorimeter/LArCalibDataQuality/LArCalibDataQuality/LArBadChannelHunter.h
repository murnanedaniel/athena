//Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


/**
 * @file  LArBadChannelHunter.h
 * @author Walter Lampl <walter.lampl @cern.ch>
 * @date Apr 2008
 * @brief Algorithm to find funny channels in LAr
 */


#ifndef LARBADCHANNELHUNTER_H
#define LARBADCHANNELHUNTER_H
 
#include <vector>
#include <string>

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "CaloIdentifier/CaloGain.h"
#include "CxxUtils/checker_macros.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "LArRecConditions/LArBadChannelCont.h"
#include "LArCabling/LArOnOffIdMapping.h"
#include "LArRecConditions/LArCalibLineMapping.h"

class LArOnlineID;
class CaloCell_ID;

class LArBadChannelHunter:public AthAlgorithm {
 
public:
  /** 
   * @brief Regular algorithm constructor
   */
  LArBadChannelHunter (const std::string& name, ISvcLocator* pSvcLocator);

  /** 
    * @brief Destructor
    */
  ~LArBadChannelHunter() {};

 /**
   * @brief Standard initialization method.
   */
  virtual StatusCode initialize() override;

  /**
    * @brief Standard execute method
    * This method has to be emtpy since all the job is done in finalize
    */
  virtual StatusCode execute() override {return StatusCode::SUCCESS;}

  /**
    * @brief Standard stop method
    */
  virtual StatusCode stop ATLAS_NOT_THREAD_SAFE () override;

  //void findFailedPatterns();

private:

  const LArOnlineID* m_onlineId; 
  const CaloCell_ID* m_caloId;

  SG::ReadCondHandleKey<LArOnOffIdMapping> m_cablingKey{this,"CablingKey","LArOnOffIdMap","SG Key of LArOnOffIdMapping object"};
  SG::ReadCondHandleKey<LArBadChannelCont> m_BCKey{this, "BadChanKey", "LArBadChannel", "SG bad channels key"};
  SG::ReadCondHandleKey<LArCalibLineMapping>  m_CLKey{this, "CalibLineKey", "LArCalibLineMap", "SG calib line key"};

  std::string m_pedKey, m_caliWaveKey;
  std::string m_avgTypeProp;
  std::string m_outFileName;
  std::string m_cutType;
  bool m_undoCorr;
  bool m_outOnlyNew;
  float m_recalcPer;

  //Thresholds:
  float m_lowNoiseTh[CaloGain::LARNGAIN]{};
  float m_highNoiseTh[CaloGain::LARNGAIN]{};
 
  float m_amplTh[CaloGain::LARNGAIN]{};
  float m_widTh[CaloGain::LARNGAIN]{};
  float m_distwidTh[CaloGain::LARNGAIN]{};
  float m_distampTh[CaloGain::LARNGAIN]{};
  float m_tmaxampTh[CaloGain::LARNGAIN]{};
  enum AvgType {
    FEB,
    PHI
  };

  AvgType m_avgType;


  const std::string channelDescription(const HWIdentifier& chid, const LArOnOffIdMapping *cabling, const unsigned gain=99) const;

  unsigned getSymId(const HWIdentifier, const LArOnOffIdMapping *cabling) const;

  class Average {
  public:
    float m_avPedRMS[CaloGain::LARNGAIN]{}; 
    float m_avPedRMSSD[CaloGain::LARNGAIN]{}; 
    unsigned m_nPed[CaloGain::LARNGAIN]{}; 
    float m_avAmpl[CaloGain::LARNGAIN]{};
    float m_avAmplSD[CaloGain::LARNGAIN]{};
    unsigned m_nAmpls[CaloGain::LARNGAIN]{};
    float m_avWid[CaloGain::LARNGAIN]{};
    float m_avWidSD[CaloGain::LARNGAIN]{};
    unsigned m_nWids[CaloGain::LARNGAIN]{};
    float m_avTmax[CaloGain::LARNGAIN]{};
    float m_avTmaxSD[CaloGain::LARNGAIN]{};
    unsigned m_nTmaxs[CaloGain::LARNGAIN]{};

    std::vector<float> m_vmedTmax[CaloGain::LARNGAIN];
    std::vector<float> m_vmedWid[CaloGain::LARNGAIN];
    std::vector<float> m_vmedAmpl[CaloGain::LARNGAIN];
    std::vector<float> m_vmedPedRMS[CaloGain::LARNGAIN];
    float m_medTmax[CaloGain::LARNGAIN]{};	
    float m_medWid[CaloGain::LARNGAIN]{};
    float m_medAmpl[CaloGain::LARNGAIN]{};
    float m_medPedRMS[CaloGain::LARNGAIN]{};
    Average();
    void finish(float); 
  };

  class CellData {
  public:
    HWIdentifier m_chid;
    float m_ampl;
    float m_wid;
    float m_tmax;	
    float m_pedRMS[3];
    CellData() : m_ampl(-1.) {
      m_wid	 =-1.;
      m_tmax	 =-1.; 		
      m_pedRMS[0]=-1.;
      m_pedRMS[1]=-1.;
      m_pedRMS[2]=-1.;
    }
  };
};



#endif

