/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "TRTDigCondBase.h"
#include "TRTDigSettings.h"
#include "TRTDigiHelper.h"

#include <cmath>
#include <cstdlib>
#include <utility>

#include "TRT_ReadoutGeometry/TRT_DetElementCollection.h"
#include "TRT_ReadoutGeometry/TRT_DetectorManager.h"
#include "TRT_ReadoutGeometry/TRT_BaseElement.h"

//helpers for identifiers and hitids:
#include "InDetIdentifier/TRT_ID.h"
#include "InDetSimEvent/TRTHitIdHelper.h"

// For the random numbers.
#include "CLHEP/Random/RandFlat.h"

//________________________________________________________________________________
TRTDigCondBase::TRTDigCondBase( const TRTDigSettings* digset,
                                const InDetDD::TRT_DetectorManager* detmgr,
                                const TRT_ID* trt_id,
                                int UseGasMix,
                                ToolHandle<ITRT_StrawStatusSummaryTool> sumTool
                                )
  : AthMessaging("TRTDigCondBase"),
    m_settings(digset),
    m_detmgr(detmgr),
    m_id_helper(trt_id),
    m_averageNoiseLevel(-1.0),
    m_crosstalk_noiselevel(-1.0),
    m_crosstalk_noiselevel_other_end(-1.0),
    m_UseGasMix(UseGasMix),
    m_sumTool(std::move(sumTool))
{
  m_crosstalk_noiselevel = m_settings->crossTalkNoiseLevel();
  m_crosstalk_noiselevel_other_end = m_settings->crossTalkNoiseLevelOtherEnd();
}


//________________________________________________________________________________
float TRTDigCondBase::strawAverageNoiseLevel() const {

  if (m_averageNoiseLevel>=0.0) {
    return m_averageNoiseLevel;
  } else {
    std::map<int,StrawState>::const_iterator it(m_hitid_to_StrawState.begin());
    m_averageNoiseLevel = 0.;
    double tmp(0.);
    for ( ; it!=m_it_hitid_to_StrawState_End; ++it ) {
      tmp += it->second.noiselevel;
    };
    m_averageNoiseLevel = tmp / totalNumberOfActiveStraws();
    return m_averageNoiseLevel;
  };
}

//________________________________________________________________________________
void TRTDigCondBase::initialize(CLHEP::HepRandomEngine* rndmEngine) {

  ATH_MSG_INFO ( "TRTDigCondBase::initialize()" );

  //id helpers:
  const TRTHitIdHelper *hitid_helper(TRTHitIdHelper::GetHelper());

  //We loop through all of the detector elements registered in the manager
  InDetDD::TRT_DetElementCollection::const_iterator it(m_detmgr->getDetectorElementBegin());
  InDetDD::TRT_DetElementCollection::const_iterator itE(m_detmgr->getDetectorElementEnd());

  unsigned int strawcount(0);
  unsigned int nBAA[3][3]  = {{0}}; // [ringwheel=0,1,2][strawGasType=0,1,2]
  unsigned int nBAC[3][3]  = {{0}}; // [ringwheel=0,1,2][strawGasType=0,1,2]
  unsigned int nECA[14][3] = {{0}}; // [ringwheel=0--13][strawGasType=0,1,2]
  unsigned int nECC[14][3] = {{0}}; // [ringwheel=0--13][strawGasType=0,1,2]

  for (;it!=itE;++it) { // loop over straws

    const double strawLength((*it)->strawLength());
    const Identifier id((*it)->identify());

    const int ringwheel(m_id_helper->layer_or_wheel(id));
    const int phisector(m_id_helper->phi_module(id));
    const int     layer(m_id_helper->straw_layer(id));
    const int      side(m_id_helper->barrel_ec(id));

    int endcap, isneg;
    switch ( side ) {
    case -2:  endcap = 1; isneg = 1; break;
    case -1:  endcap = 0; isneg = 1; break;
    case  1:  endcap = 0; isneg = 0; break;
    case  2:  endcap = 1; isneg = 0; break;
    default:
      std::cout << "TRTDigitization::TRTDigCondBase::createListOfValidStraws "
                << "FATAL - identifier problems - skipping detector element!!" << std::endl;
      continue;
    }

    for (unsigned int iStraw(0); iStraw <(*it)->nStraws(); ++iStraw) {

      // get ID of the straw, and the gas mix
      const int hitid(hitid_helper->buildHitId( endcap, isneg, ringwheel, phisector, layer, iStraw));
      Identifier strawId = m_id_helper->straw_id(side, phisector, ringwheel, layer, iStraw);
      const int statusHT = m_sumTool->getStatusHT(strawId);
      const int strawGasType = TRTDigiHelper::StrawGasType(statusHT,m_UseGasMix, &msg());

      //Get info about the straw conditions, then create and fill the strawstate
      double noiselevel, relative_noiseamplitude;
      setStrawStateInfo( strawId, strawLength, noiselevel, relative_noiseamplitude, rndmEngine );
      StrawState strawstate{};
      strawstate.noiselevel = noiselevel; // same for all gas types
      strawstate.lowthreshold = ( !(hitid & 0x00200000) ) ? m_settings->lowThresholdBar(strawGasType) : m_settings->lowThresholdEC(strawGasType);
      strawstate.noiseamplitude= relative_noiseamplitude; // These two are later regulated noise code
      m_hitid_to_StrawState[ hitid ] = strawstate; // Insert into the map:

      // Count the number of straws
      ++strawcount;

      // Count the gas fraction in a number of regions:
      if ( side==+1 && ringwheel>=0 && ringwheel<=2  ) nBAA[ringwheel][strawGasType]++; // [ringwheel=0,1,2][strawGasType=0,1,2]
      if ( side==-1 && ringwheel>=0 && ringwheel<=2  ) nBAC[ringwheel][strawGasType]++; // [ringwheel=0,1,2][strawGasType=0,1,2]
      if ( side==+2 && ringwheel>=0 && ringwheel<=13 ) nECA[ringwheel][strawGasType]++; // [ringwheel=0, 13][strawGasType=0,1,2]
      if ( side==-2 && ringwheel>=0 && ringwheel<=13 ) nECC[ringwheel][strawGasType]++; // [ringwheel=0, 13][strawGasType=0,1,2]

    };

  }; // end "loop over straws"

  m_it_hitid_to_StrawState = m_hitid_to_StrawState.begin();
  m_it_hitid_to_StrawState_End = m_hitid_to_StrawState.end();

  //just put it to something:
  m_it_hitid_to_StrawState_Last = m_it_hitid_to_StrawState;

  if (m_hitid_to_StrawState.empty()) {
    ATH_MSG_ERROR("TRTDigCondBase::initialize it seems that ALL straws are dead/masked! This wont work.");
  }

  //just to avoid having an uninitialized iterator hanging around:
  m_all_it_hitid_to_StrawState_previous = m_hitid_to_StrawState.begin();

  // Finally give some useful information about the gas mix chosen and that which we actually get!
  if ( m_UseGasMix==0) std::cout << "TRTDigCondBase  INFO Gas Geometry: UseGasMix==0; using StatusHT to determine the gas geometry." << std::endl;
  if ( m_UseGasMix==1) std::cout << "TRTDigCondBase  INFO Gas Geometry: UseGasMix==1; expect Xenon in the entire detector." << std::endl;
  if ( m_UseGasMix==2) std::cout << "TRTDigCondBase  INFO Gas Geometry: UseGasMix==2; expect Krypton in the entire detector." << std::endl;
  if ( m_UseGasMix==3) std::cout << "TRTDigCondBase  INFO Gas Geometry: UseGasMix==3; expect Argon in the entire detector." << std::endl;
  if ( m_UseGasMix<0 || m_UseGasMix>3) {
    std::cout << "TRTDigCondBase  ERROR Gas Geometry: UseGasMix==" << m_UseGasMix << ", must be 0,1,2or 3!" << std::endl;
    throw std::exception();
  }
  std::cout << "TRTDigCondBase  INFO Gas Geometry: strawcount=" << strawcount << std::endl;

  std::cout << "TRTDigCondBase  INFO Gas Geometry: BA_A[Xe] = " << nBAA[0][0] << " " << nBAA[1][0] << " " << nBAA[2][0] << std::endl;
  std::cout << "TRTDigCondBase  INFO Gas Geometry: BA_A[Kr] = " << nBAA[0][1] << " " << nBAA[1][1] << " " << nBAA[2][1] << std::endl;
  std::cout << "TRTDigCondBase  INFO Gas Geometry: BA_A[Ar] = " << nBAA[0][2] << " " << nBAA[1][2] << " " << nBAA[2][2] << std::endl;

  std::cout << "TRTDigCondBase  INFO Gas Geometry: BA_C[Xe] = " << nBAC[0][0] << " " << nBAC[1][0] << " " << nBAC[2][0] << std::endl;
  std::cout << "TRTDigCondBase  INFO Gas Geometry: BA_C[Kr] = " << nBAC[0][1] << " " << nBAC[1][1] << " " << nBAC[2][1] << std::endl;
  std::cout << "TRTDigCondBase  INFO Gas Geometry: BA_C[Ar] = " << nBAC[0][2] << " " << nBAC[1][2] << " " << nBAC[2][2] << std::endl;

  std::cout << "TRTDigCondBase  INFO Gas Geometry: EC_A[Xe] = "; for (int i=0;i<14;i++) std::cout << nECA[i][0] << " "; std::cout << std::endl;
  std::cout << "TRTDigCondBase  INFO Gas Geometry: EC_A[Kr] = "; for (int i=0;i<14;i++) std::cout << nECA[i][1] << " "; std::cout << std::endl;
  std::cout << "TRTDigCondBase  INFO Gas Geometry: EC_A[Ar] = "; for (int i=0;i<14;i++) std::cout << nECA[i][2] << " "; std::cout << std::endl;

  std::cout << "TRTDigCondBase  INFO Gas Geometry: EC_C[Xe] = "; for (int i=0;i<14;i++) std::cout << nECC[i][0] << " "; std::cout << std::endl;
  std::cout << "TRTDigCondBase  INFO Gas Geometry: EC_C[Kr] = "; for (int i=0;i<14;i++) std::cout << nECC[i][1] << " "; std::cout << std::endl;
  std::cout << "TRTDigCondBase  INFO Gas Geometry: EC_C[Ar] = "; for (int i=0;i<14;i++) std::cout << nECC[i][2] << " "; std::cout << std::endl;

}

//________________________________________________________________________________
void TRTDigCondBase::resetGetNextNoisyStraw() {
  m_it_hitid_to_StrawState = m_hitid_to_StrawState.begin();
}

//________________________________________________________________________________
bool TRTDigCondBase::getNextNoisyStraw( CLHEP::HepRandomEngine* randengine, int& hitID, float& noiselvl ) {
  for (;m_it_hitid_to_StrawState!=m_it_hitid_to_StrawState_End;++m_it_hitid_to_StrawState) {
    noiselvl = m_it_hitid_to_StrawState->second.noiselevel;
    if ( CLHEP::RandFlat::shoot(randengine, 0.0, 1.0) < noiselvl ) {
      ++m_it_hitid_to_StrawState; //Important! if removed, iterator is not incremented in case rand<noiselevel!!!
      hitID = m_it_hitid_to_StrawState->first;
      return true;
    };
  };
  return false;
}

//________________________________________________________________________________
bool TRTDigCondBase::crossTalkNoise( CLHEP::HepRandomEngine* randengine ) {
  const float noise(- m_crosstalk_noiselevel * log(CLHEP::RandFlat::shoot(randengine, 0.0, 1.0)));
  if ( CLHEP::RandFlat::shoot(randengine, 0.0, 1.0) < noise ) return true;
  return false;
}

//________________________________________________________________________________
bool TRTDigCondBase::crossTalkNoiseOtherEnd( CLHEP::HepRandomEngine* randengine ) {
  const float noise(- m_crosstalk_noiselevel_other_end * log(CLHEP::RandFlat::shoot(randengine, 0.0, 1.0)));
  if ( CLHEP::RandFlat::shoot(randengine, 0.0, 1.0) < noise ) return true;
  return false;
}
