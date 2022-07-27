/*
Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
#ifndef TRIGT2DJTRIG_ED_HYPOTOOL_H
#define TRIGT2DJTRIG_ED_HYPOTOOL_H

#include "Gaudi/Property.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "TrigCompositeUtils/HLTIdentifier.h"
#include "TrigCompositeUtils/TrigCompositeUtils.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODJet/Jet.h"
#include "AthContainers/AuxElement.h"
#include <vector>
#include "AthenaMonitoringKernel/Monitored.h"

class DisplacedJetEventDecisionHypoTool : virtual public ::AthAlgTool
{
public:

  DisplacedJetEventDecisionHypoTool( const std::string& type,
    const std::string& name,
    const IInterface* parent );

    virtual ~DisplacedJetEventDecisionHypoTool() = default;
    virtual StatusCode initialize() override;

    struct DecisionTuple {
      TrigCompositeUtils::Decision* outDecision;
      const TrigCompositeUtils::Decision* inDecision;
      const TrigCompositeUtils::DecisionIDContainer previousDecisionIDs;
      const xAOD::Jet* jet;
      const xAOD::TrigComposite* counts;
      const xAOD::TrigComposite* info;
    };

    StatusCode decide( std::vector<DecisionTuple>& decs ) const;

  private:

    HLT::Identifier m_decisionId;

    Gaudi::Property<int> m_min_h_jets{this, "min_h_jets",{1}, "m_min_h_jets"};
    Gaudi::Property<int> m_min_l_jets{this, "min_l_jets",{1}, "m_min_h_jets"};
    
    Gaudi::Property<std::string> m_cutname{this, "cut_name",{""}, "Name of cuts, used for decoration names"};

    ToolHandle<GenericMonitoringTool> m_monTool{this,"MonTool","","Monitoring tool"};
  };

  #endif //> !TRIGT2MINBIAS_TRACKCOUNTHYPOTOOL_H
