
/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGBJETHYPO_TRIGBJETBTAGHYPOALG_H
#define TRIGBJETHYPO_TRIGBJETBTAGHYPOALG_H 1

#include "TrigBjetHypoAlgBase.h"
#include "TrigBjetBtagHypoTool.h"

#include <string>

#include "TrigCompositeUtils/TrigCompositeUtils.h"

#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"

#include "xAODBTagging/BTaggingAuxContainer.h"
#include "xAODBTagging/BTaggingContainer.h"

#include "AthenaMonitoringKernel/Monitored.h"
#include "AthenaMonitoringKernel/GenericMonitoringTool.h"

#define MONITOR_BTAG_AUX_TRACK_VAR(VAR_NAME, VAR_TYPE ) \
      auto monitor_for_##VAR_NAME = Monitored::Collection( #VAR_NAME, \
        (*bTagLink)->auxdata< std::vector<VAR_TYPE> >( #VAR_NAME ));

#define MONITOR_BTAG_AUX_VAR(VAR_NAME, VAR_TYPE, CONTAINER ) \
    auto monitor_for_##VAR_NAME = Monitored::Collection( #VAR_NAME, CONTAINER, \
      [](const ElementLink< xAOD::BTaggingContainer >& bTagLink) {  \
        return (*bTagLink)->auxdata<VAR_TYPE>( #VAR_NAME ); \
      } \
    );


class TrigBjetBtagHypoAlg : public TrigBjetHypoAlgBase {
 public:
  TrigBjetBtagHypoAlg( const std::string& name, ISvcLocator* pSvcLocator );

  virtual StatusCode  initialize();
  virtual StatusCode  execute( const EventContext& context ) const;

 private:
  TrigBjetBtagHypoAlg();

  // online monitoring 
  virtual StatusCode monitor_jets( const ElementLinkVector<xAOD::JetContainer >& jetELs, const ElementLinkVector<xAOD::JetContainer >& all_bTaggedJetELs ) const ;
  virtual StatusCode monitor_tracks( const EventContext& context, const TrigCompositeUtils::DecisionContainer* prevDecisionContainer ) const;
  virtual StatusCode monitor_primary_vertex( const ElementLink< xAOD::VertexContainer >& primVertexEL ) const;
  virtual StatusCode monitor_flavor_probabilities( const ElementLinkVector< xAOD::BTaggingContainer >& bTaggingEL, const std::string& var_name) const;
  virtual StatusCode monitor_flavor_bb_probabilities( const ElementLinkVector< xAOD::BTaggingContainer >& bTaggingEL, const std::string& var_name) const;
  virtual ElementLinkVector<xAOD::BTaggingContainer> collect_valid_links(
      const ElementLinkVector< xAOD::BTaggingContainer >& bTaggingEL, std::string tagger ) const;
  virtual StatusCode monitor_btagging( const ElementLinkVector< xAOD::BTaggingContainer >& bTaggingEL ) const;
  
 private:
  ToolHandleArray< TrigBjetBtagHypoTool > m_hypoTools {this,"HypoTools",{},"Hypo Tools"};
  ToolHandle<GenericMonitoringTool> m_monTool{this,"MonTool","","Monitoring tool"};
  
  SG::ReadHandleKey< xAOD::JetContainer > m_bTaggedJetKey {this,"BTaggedJetKey","","Key for b-tagged jets"};
  SG::ReadHandleKey< xAOD::BTaggingContainer> m_bTagKey {this,"BTaggingKey","","Key for BTagging"};
  SG::ReadHandleKey< xAOD::TrackParticleContainer > m_trackKey {this,"TracksKey","","Key for precision tracks"};
  SG::ReadHandleKey< xAOD::VertexContainer > m_inputPrmVtx {this,"PrmVtxKey","","Key for Primary vertex collection for monitoring"};

  Gaudi::Property< std::string > m_bTaggingLink {this,"BTaggingLink","Unspecified","b-Tagging Link name in navigation (output)"};
  Gaudi::Property< std::string > m_prmVtxLink {this,"PrmVtxLink","Unspecified","Vertex Link name in navigation (input)"};
  Gaudi::Property<std::string> m_btaggingLinkName{this, "BtaggingLinkName", "btag"}; // TM 2021-10-30

  SG::ReadCondHandleKey< InDet::BeamSpotData > m_beamSpotKey{ this,
     "BeamSpotKey", "BeamSpotData", "SG key for beam spot" };

};

#endif

