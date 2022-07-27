/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
#ifndef TRIGOUTPUTHANDLING_STREAMTAGMAKERTOOL_H
#define TRIGOUTPUTHANDLING_STREAMTAGMAKERTOOL_H

// Trigger includes
#include "TrigConfData/HLTMenu.h"
#include "TrigCompositeUtils/TrigCompositeUtils.h"
#include "TrigOutputHandling/HLTResultMTMakerTool.h"
#include "TrigPartialEventBuilding/PEBInfoWriterToolBase.h" // Defines PEBInfo and keys to retrieve it

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"

// System includes
#include <string>
#include <unordered_map>
#include <tuple>

/**
 * @class StreamTagMakerTool
 * @brief Makes stream tags out of chain decisions
 **/
class StreamTagMakerTool : public extends<AthAlgTool, HLTResultMTMakerTool> {
public:
  StreamTagMakerTool(const std::string& type, const std::string& name, const IInterface* parent);
  virtual ~StreamTagMakerTool() override = default;

  virtual StatusCode fill( HLT::HLTResultMT& resultToFill, const EventContext& ctx ) const override;

  virtual StatusCode initialize() override;
  virtual StatusCode start() override;

  /// Type describing StreamTag information needed by the tool: {name, type, obeysLumiBlock, forceFullEventBuilding}
  using StreamTagInfo = std::tuple<std::string, std::string, bool, bool>;

private:
  SG::ReadHandleKey<TrigConf::HLTMenu> m_hltMenuKey{"DetectorStore+HLTTriggerMenu"};

  SG::ReadHandleKey<TrigCompositeUtils::DecisionContainer> m_finalChainDecisions {this, "ChainDecisions", "HLTNav_Summary",
    "Container with final chain decisions" };

  SG::ReadHandleKeyArray<TrigCompositeUtils::DecisionContainer> m_pebDecisionKeys {this, "PEBDecisionKeys", {},
    "Decisions including PEBInfo" };

  /// Chain to streams map filled from the HLT Menu JSON
  std::unordered_map<TrigCompositeUtils::DecisionID, std::vector<StreamTagInfo> > m_mapping;

  /// Helper method to fill the chainID->PEBInfo map
  StatusCode fillPEBInfoMap(std::unordered_map<TrigCompositeUtils::DecisionID, PEBInfoWriterToolBase::PEBInfo>& map, const EventContext& ctx) const;
};


#endif // TRIGOUTPUTHANDLING_STREAMTAGMAKERTOOL_H
