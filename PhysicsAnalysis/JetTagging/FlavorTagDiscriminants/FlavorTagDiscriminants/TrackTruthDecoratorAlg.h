/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRACK_TRUTH_DECORATOR_ALG_HH
#define TRACK_TRUTH_DECORATOR_ALG_HH

#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "AthContainers/AuxElement.h"
#include "StoreGate/WriteDecorHandleKey.h"
#include "StoreGate/ReadDecorHandleKey.h"

#include "xAODTracking/TrackParticleContainerFwd.h"
#include "InDetTrackSystematicsTools/InDetTrackTruthOriginTool.h"
#include "xAODTruth/TruthEventContainer.h"


namespace FlavorTagDiscriminants {

  class TrackTruthDecoratorAlg: public AthReentrantAlgorithm {
  public:
    TrackTruthDecoratorAlg(const std::string& name,
                          ISvcLocator* pSvcLocator );

    virtual StatusCode initialize() override;
    virtual StatusCode execute(const EventContext& ) const override;

  private:

    // Input Containers
    SG::ReadHandleKey< xAOD::TrackParticleContainer > m_TrackContainerKey {
      this,"trackContainer", "InDetTrackParticles",
        "Key for the input track collection"};
    SG::ReadHandleKey< xAOD::TruthEventContainer > m_TruthEventsKey {
      this,"truthEvents", "TruthEvents",
        "Key for the input truth event collection"};

    // Decorators for tracks
    SG::WriteDecorHandleKey< xAOD::TrackParticleContainer > m_dec_origin_label {
      this, "truthOriginLabel", "truthOriginLabel", 
        "Exclusive origin label of the track"};
    SG::WriteDecorHandleKey< xAOD::TrackParticleContainer > m_dec_type_label {
      this, "truthTypeLabel", "truthTypeLabel", 
        "Exclusive truth type label of the track"};
    SG::WriteDecorHandleKey< xAOD::TrackParticleContainer > m_dec_vertex_index {
      this, "truthVertexIndex", "truthVertexIndex", 
        "Truth vertex index of the track"};
    SG::WriteDecorHandleKey< xAOD::TrackParticleContainer > m_dec_barcode {
      this, "barcode", "barcode", 
        "Barcode of linked truth particle"};
    SG::WriteDecorHandleKey< xAOD::TrackParticleContainer > m_dec_parent_barcode {
      this, "parentBarcode", "parentBarcode", 
        "Barcode of parent of linked truth particle"};

    // truth origin tool
    ToolHandle<InDet::InDetTrackTruthOriginTool> m_trackTruthOriginTool {
      this, "trackTruthOriginTool", "InDet::InDetTrackTruthOriginTool", 
        "track truth origin tool"};

    Gaudi::Property<float> m_truthVertexMergeDistance {
      this, "truthVertexMergeDistance", 0.1, 
        "Merge any truth vertices within this distance [mm]"};
  };

}

#endif
