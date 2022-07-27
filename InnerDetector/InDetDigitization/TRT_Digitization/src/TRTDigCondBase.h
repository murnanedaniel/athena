/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRT_DIGITIZATION_TRTDIGCONDBASE_H
#define TRT_DIGITIZATION_TRTDIGCONDBASE_H

#include "AthenaBaseComps/AthMessaging.h"
#include "CxxUtils/checker_macros.h"
#include "Identifier/Identifier.h"
#include "TRT_ConditionsServices/ITRT_StrawStatusSummaryTool.h" // added by Sasha for Argon

#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "CLHEP/Random/RandomEngine.h"

#include <atomic>
#include <map>
#include <mutex>
#include <set>

class TRT_ID;

namespace InDetDD {
  class TRT_DetectorManager;
}
class TRTDigSettings;
/**
 * Communication with CondDB
 */
class TRTDigCondBase : public AthMessaging {

public:

  /** Constructor */
  TRTDigCondBase( const TRTDigSettings*,
                  const InDetDD::TRT_DetectorManager*,
                  const TRT_ID*,
                  int UseGasMix,
                  ToolHandle<ITRT_StrawStatusSummaryTool> sumTool
                  );

  /** Destructor */
  virtual ~TRTDigCondBase() {};

  /**
   * @note Must be called exactly once before the class is used for anything.
   */
  void initialize(CLHEP::HepRandomEngine* rndmEngine);

  /** Get average noise level in straw */
  float strawAverageNoiseLevel() const;         //[first time slow, else fast]

  /** Get total number of active straws */
  unsigned int totalNumberOfActiveStraws() const; //[fast]

  /**
   * Get straw data
   * mixed condition is implemented
   * function will return both Argon and Xenon straws (with according LT)
   * for Argon LT ~300eV, for Argon ~100eV
   * @param hitID:          straw ID
   * @param lowthreshold:   low threshold discriminator level
   * @param noiseamplitude: noise amplitude
   */
  void getStrawData( const int& hitID,
                     double& lowthreshold,
                     double& noiseamplitude ) const;

  //--- For noise in unhit straws:

  /**
   * For noise in unhit straws: Rewind straw list to start from beginning
   */
  void resetGetNextNoisyStraw();

  /**
   * For simulation of noise in unhit straws: get next noisy straw.
   * @param hitID: ID of next noisy straw
   * @param noiselevel: noise level
   */
  bool getNextNoisyStraw( CLHEP::HepRandomEngine*, int& hitID, float& noiselevel );


  //Crosstalk noise
  bool crossTalkNoise( CLHEP::HepRandomEngine* );
  bool crossTalkNoiseOtherEnd( CLHEP::HepRandomEngine* );

  //--- For looping over all straws (only to be used at initialization):

  /**
   * For looping over all straws: Rewind straw list to start from beginning
   **/
  void resetGetNextStraw();

  /**
   * For looping over all straws: get next straw
   * @param hitID: ID of next straw
   * @param lowthreshold: low threshold discrimenator setting
   * @param noiseamplitude: noise amplitude
   */
  bool getNextStraw( int& hitID, float& noiselevel, float& noiseamp );

  /**
   * Set straw parameters
   * @param hitID:          ID of straw
   * @param lowthreshold:   low threshold discriminator setting
   * @param noiseamplitude: noise amplitude
   */
  void setRefinedStrawParameters( const int& hitID,
                                  const double& lowthreshold,
                                  const double& noiseamplitude );

protected:

  /**
   * Get straw state info based on hitid and strawlength.
   *
   * @note  This must be overriden by the actual map provider
   * @param TRT_Identifier: Identifier of straw (input)
   * @param strawlength: straw length (input)
   * @param noiselevel: noise level this straw (e.g. 0.02 for 2%)
   * @param relative_noiseamplitude: the relative noise amplitude as compared
   *                                 to other straws. The overall normalization
   *                                 doesn't matter as all noise amplitudes
   *                                 will be rescaled by the same factor.
   *                                 Should never be put to zero!
   */
  virtual void setStrawStateInfo(Identifier& TRT_Identifier,
                                 const double& strawlength,
                                 double& noiselevel,
                                 double& relative_noiseamplitude,
                                 CLHEP::HepRandomEngine *rndmEngine) = 0;

protected:
  const TRTDigSettings* m_settings;
  const InDetDD::TRT_DetectorManager* m_detmgr;
  const TRT_ID* m_id_helper;

private:
  /** Straw state */
  struct StrawState {
    float noiselevel;     /**< Noise level                         */
    float noiseamplitude; /**< Noise amplitude                     */
    float lowthreshold;   /**< Low threshold discriminator setting */
  };

  //--- Data maps:

  /** Global map from straw ID to straw state */
  // StrawState is struct with members:
  // 1) noiselevel      - the same for Xenon/Argon straws
  // 2) noiseamplitude  - the same for Xenon/Argon straws
  // 3) lowthreshold    - different for Xenon and Argon straws
  // for more info see info for setStrawStateInfo function, 28 lines above!

  std::map<int,StrawState> m_hitid_to_StrawState;  //<--------- 7MB!!!!

  //--- iterators:

  /** Iterator over straw state map */
  std::map<int,StrawState>::const_iterator m_it_hitid_to_StrawState;
  /** Iterator pointing to last straw in straw state map */
  std::map<int,StrawState>::const_iterator m_it_hitid_to_StrawState_End;
  /** Iterator used for caching (ought to be called _Previous) */
  mutable std::map<int,StrawState>::const_iterator m_it_hitid_to_StrawState_Last ATLAS_THREAD_SAFE; // Guarded by m_mutex
  mutable std::mutex m_mutex;

  mutable std::atomic<float> m_averageNoiseLevel; /**< Average noise level */
  double m_crosstalk_noiselevel;
  double m_crosstalk_noiselevel_other_end;

  //---  For iterating through all of the straw hitids:
  /** Iterator over straw state map */
  std::map<int,StrawState>::iterator m_all_it_hitid_to_StrawState;
  /** Iterator used for caching */
  std::map<int,StrawState>::iterator m_all_it_hitid_to_StrawState_previous;

protected:
  int m_UseGasMix;
  ToolHandle<ITRT_StrawStatusSummaryTool> m_sumTool; // added by Sasha for Argon

};

////////////////
/// INLINES ////
////////////////

//___________________________________________________________________________
inline void TRTDigCondBase::getStrawData( const int& hitID,
                                          double& lowthreshold,
                                          double& noiseamplitude ) const {
  //fixme: do we actually benefit from this caching?
  std::lock_guard<std::mutex> lock(m_mutex);
  if (m_it_hitid_to_StrawState_Last->first != hitID)
    m_it_hitid_to_StrawState_Last = m_hitid_to_StrawState.find(hitID);

  lowthreshold   = m_it_hitid_to_StrawState_Last->second.lowthreshold;
  noiseamplitude = m_it_hitid_to_StrawState_Last->second.noiseamplitude;
}

//___________________________________________________________________________
inline unsigned int TRTDigCondBase::totalNumberOfActiveStraws() const {
  return m_hitid_to_StrawState.size();
}

//___________________________________________________________________________
inline void TRTDigCondBase::resetGetNextStraw() {
  m_all_it_hitid_to_StrawState = m_hitid_to_StrawState.begin() ;
}

//___________________________________________________________________________
inline bool TRTDigCondBase::getNextStraw( int& hitID,
                                          float& noiselevel,
                                          float& noiseamp ) {
  m_all_it_hitid_to_StrawState_previous = m_all_it_hitid_to_StrawState;
  if ( ++m_all_it_hitid_to_StrawState == m_it_hitid_to_StrawState_End) {
    hitID      = 0;
    noiselevel = 0.;
    return false;
  } else {
    hitID      = m_all_it_hitid_to_StrawState->first;
    noiselevel = m_all_it_hitid_to_StrawState->second.noiselevel;
    noiseamp   = m_all_it_hitid_to_StrawState->second.noiseamplitude;
    return true;
  };
}

//___________________________________________________________________________
inline void TRTDigCondBase::setRefinedStrawParameters( const int& hitID,
                                                       const double& lowthreshold,
                                                       const double& noiseamplitude ) {
  //we might have the iterator cached already (fixme... nah)
  if (hitID == m_all_it_hitid_to_StrawState_previous->first) {
    m_all_it_hitid_to_StrawState_previous->second.lowthreshold = lowthreshold;
    m_all_it_hitid_to_StrawState_previous->second.noiseamplitude =
      noiseamplitude;
  } else {
    m_all_it_hitid_to_StrawState_previous = m_hitid_to_StrawState.find(hitID);
    m_hitid_to_StrawState.find(hitID)->second.lowthreshold = lowthreshold;
    m_hitid_to_StrawState.find(hitID)->second.noiseamplitude = noiseamplitude;
  };
}

#endif
