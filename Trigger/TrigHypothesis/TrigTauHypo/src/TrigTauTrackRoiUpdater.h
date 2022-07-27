/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGTAUHYPO_TrigTauTrackRoiUpdater_H
#define TRIGTAUHYPO_TrigTauTrackRoiUpdater_H

#include <iostream>

#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "TrkTrack/TrackCollection.h"

#include "TrigSteeringEvent/TrigRoiDescriptor.h"

#include "xAODTau/TauJetContainer.h"
#include "tauRecTools/BDTHelper.h"

class TrigTauTrackRoiUpdater : public AthReentrantAlgorithm {

 public:
  TrigTauTrackRoiUpdater(const std::string&, ISvcLocator*);
  ~TrigTauTrackRoiUpdater();
  
  virtual StatusCode initialize() override;
  virtual StatusCode execute(const EventContext&) const override;

 private:

  Gaudi::Property< float > m_z0HalfWidth  {this,"z0HalfWidth",7.0,"z0 Half width for FTF Iso"};
  Gaudi::Property< float > m_etaHalfWidth {this,"etaHalfWidth",0.4,"eta Half width for FTF Iso"};
  Gaudi::Property< float > m_phiHalfWidth {this,"phiHalfWidth",0.4,"phi Half width for FTF Iso"};
  Gaudi::Property< int > m_nHitPix {this,"nHitPix",2,"at least n hits in pixels on lead track"};
  Gaudi::Property< int > m_nSiHoles {this,"nSiHoles",2,"maximum number of Si holes on lead track"};
  Gaudi::Property< std::string > m_BDTweights {this,"BDTweights","","path to BDT file, when tauIso ROI is defined by highest-BDT-score tauCore track (empty string means BDT is not used)"};

  std::unique_ptr<tauRecTools::BDTHelper> m_reader;

  struct BDTInputVariables
  {
    float logtrk_pt{0.0};
    float abstrck_z0{0.0};
    float abstrk_d0{0.0};
    float trk_nPiHits{0.0};
    float trk_nSiHoles{0.0};
    float logtrk_ratiopt{0.0};
    float trk_dR{0.0};
    float trk_dRtoleadtrk{0.0};
    float CaloHad_pt{0.0};
    float CaloEM_pt{0.0};
  };

  SG::ReadHandleKey< TrigRoiDescriptorCollection > m_roIInputKey {this,"RoIInputKey","InputRoI","Input RoI key name"};
  SG::ReadHandleKey< TrackCollection > m_tracksKey { this, "fastTracksKey", "fasttracks", "fast tracks in view" };
  SG::WriteHandleKey< TrigRoiDescriptorCollection > m_roIOutputKey {this,"RoIOutputKey","InViewRoI","Output RoI Collection Key"};

  SG::ReadHandleKey< xAOD::TauJetContainer> m_tauJetKey      { this, "Key_trigTauJetInputContainer", "HLT_taujet", "input taujet container" };
  double getBDTscore(const xAOD::TauJet* tau, const Trk::Track* track, const Trk::Track* leadtrack ) const;

};
#endif
