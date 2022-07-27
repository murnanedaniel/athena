/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "StreamTagMakerTool.h"
#include "TrigConfData/HLTMenu.h"
#include "TrigConfData/HLTChain.h"
#include "GaudiKernel/IAlgExecStateSvc.h"
#include "eformat/StreamTag.h"
#include <sstream>

using namespace TrigCompositeUtils;

namespace {
  std::string formatStreamTagInfo(const StreamTagMakerTool::StreamTagInfo& info) {
    std::ostringstream ss;
    ss << "[" << std::get<0>(info) << ", " << std::get<1>(info) << ", " << std::get<2>(info) << ", " << std::get<3>(info) << "]";
    return ss.str();
  }
}

// =============================================================================

StreamTagMakerTool::StreamTagMakerTool( const std::string& type, const std::string& name, const IInterface* parent ):
  base_class( type, name, parent ) {}

// =============================================================================

StatusCode StreamTagMakerTool::initialize() {
  renounceArray( m_pebDecisionKeys );
  ATH_CHECK( m_pebDecisionKeys.initialize() );
  ATH_CHECK( m_finalChainDecisions.initialize() );
  ATH_CHECK( m_hltMenuKey.initialize() );
  return StatusCode::SUCCESS;
}

// =============================================================================

StatusCode StreamTagMakerTool::start() {
  m_mapping.clear();
  auto hltMenu = SG::makeHandle( m_hltMenuKey );
  if( ! hltMenu.isValid() ) {
    ATH_MSG_FATAL("Failed to get the HLT menu from the DetectorStore");
    return StatusCode::FAILURE;
  }

  ATH_MSG_INFO("Configuring from HLTMenu from DetStore with " << hltMenu->size() << " chains");

  std::vector<TrigConf::DataStructure> allStreams = hltMenu->streams();
  ATH_MSG_INFO("Menu has " << allStreams.size() << " streams defined");
  std::map<std::string, StreamTagInfo> streamDictionary;
  ATH_MSG_DEBUG("StreamTags [Name,type,obeyLB,forceFullEventBuilding]:");
  for (const TrigConf::DataStructure& stream : allStreams) {
    try {
      std::string stream_name         = stream.getAttribute("name");
      std::string stream_type         = stream.getAttribute("type");
      std::string_view obeyLB         = stream.getAttribute("obeyLB");
      std::string_view fullEventBuild = stream.getAttribute("forceFullEventBuilding");
      StreamTagInfo streamTag = {
        stream_name,
        stream_type,
        obeyLB == "true",
        fullEventBuild == "true"
      };
      ATH_MSG_DEBUG("-- " << formatStreamTagInfo(streamTag));
      streamDictionary.insert(std::pair<std::string, StreamTagInfo>(stream.getAttribute("name"),streamTag));
    } catch (const std::exception& ex) {
      ATH_MSG_ERROR("Failure reading stream tag configuration from JSON: " << ex.what());
      return StatusCode::FAILURE;
    }
  }

  for (const TrigConf::Chain & chain : *hltMenu) {
    std::vector<std::string> streams = chain.streams();
    if (streams.empty()) {
      ATH_MSG_ERROR("Chain " << chain.name() << " has no streams assigned");
      return StatusCode::FAILURE;
    }
    ATH_MSG_DEBUG("Chain " << chain.name() << " is assigned to " << streams.size() << " streams");
    m_mapping[ HLT::Identifier( chain.name() ) ] = {};
    for (const std::string& stream : streams) {
      ATH_MSG_DEBUG("-- " << stream);
      if (const auto streamIt = streamDictionary.find(stream); streamIt != streamDictionary.end()) {
         m_mapping[ HLT::Identifier(chain.name()).numeric() ].push_back(streamIt->second);
      }else{
         ATH_MSG_ERROR("Failure reading stream tag configuration for stream: " << stream);
         return StatusCode::FAILURE;
      }
    }
  }
  return StatusCode::SUCCESS;
}

// =============================================================================

StatusCode StreamTagMakerTool::fill( HLT::HLTResultMT& resultToFill, const EventContext& ctx ) const {
  // obtain chain decisions,
  using namespace TrigCompositeUtils;
  SG::ReadHandle<DecisionContainer> chainsHandle{m_finalChainDecisions, ctx};
  if (!chainsHandle.isValid()) {
    SmartIF<IAlgExecStateSvc> aess = svcLoc()->service<IAlgExecStateSvc>("AlgExecStateSvc", false);
    if (aess.isValid() && aess->eventStatus(ctx) != EventStatus::Success) {
      ATH_MSG_WARNING("Failed event, " << m_finalChainDecisions.key() << " is unavailable. Skipping stream tag making.");
      return StatusCode::SUCCESS;
    }
    ATH_MSG_ERROR("Unable to read in the " << m_finalChainDecisions.key() << " from the DecisionSummaryMakerAlg");
    return StatusCode::FAILURE;
  }

  // Extract DecisionIDContainer from DecisionContainer
  auto getPassIDs = [this](const DecisionContainer& finalDecisions, const std::string& nodeName, DecisionIDContainer& passIDs) {
    const Decision* passNode = getNodeByName(finalDecisions, nodeName);
    if (passNode==nullptr) {
      ATH_MSG_ERROR("Unable to read in the " << nodeName << " node from the " << m_finalChainDecisions.key() << " collection from the DecisionSummaryMakerAlg");
      return StatusCode::FAILURE;
    }
    decisionIDs(passNode, passIDs);
    return StatusCode::SUCCESS;
  };

  // Get raw pass decisions
  DecisionIDContainer passRawIDs;
  ATH_CHECK(getPassIDs(*chainsHandle,summaryPassNodeName(),passRawIDs));
  if (passRawIDs.empty()) {
    ATH_MSG_DEBUG("No chains passed, event rejected");
    return StatusCode::SUCCESS;
  }

  // Get express pass decisions
  DecisionIDContainer passExpressIDs;
  ATH_CHECK(getPassIDs(*chainsHandle,summaryPassExpressNodeName(),passExpressIDs));

  std::unordered_map<unsigned int, PEBInfoWriterToolBase::PEBInfo> chainToPEBInfo;
  ATH_CHECK(fillPEBInfoMap(chainToPEBInfo, ctx));

  // for each accepted chain look up the map of chainID -> ST
  for ( DecisionID chain: passRawIDs ) {

    if (TrigCompositeUtils::isLegId(HLT::Identifier(chain)) )
      continue;
    
    auto mappingIter = m_mapping.find( chain );
    if( mappingIter == m_mapping.end() ) {
      ATH_MSG_ERROR("Each chain has to have the StreamTag associated whereas the " << HLT::Identifier( chain ) << " does not" );
      return StatusCode::FAILURE;
    }
    
    const std::vector<StreamTagInfo>& streams = mappingIter->second;
    for (const StreamTagInfo& streamTagInfo : streams) {
      const auto& [st_name, st_type, obeysLB, forceFullEvent] = streamTagInfo;
      ATH_MSG_DEBUG("Chain " << HLT::Identifier( chain ) << " accepted event into stream " << st_type << "_" << st_name
                    << " (obeysLB=" << obeysLB << ", forceFullEvent=" << forceFullEvent << ")");

      // express stream prescaling
      if (st_type == "express") {
        const auto activeExpress = std::find(passExpressIDs.begin(), passExpressIDs.end(), chain);
        // chain didn't pass express stream prescale so don't add stream tag
        if (activeExpress == passExpressIDs.end()) continue;
      }

      std::set<uint32_t> robs;
      std::set<eformat::SubDetector> subdets;
      if (!forceFullEvent) {
        auto it = chainToPEBInfo.find(chain);
        if (it != chainToPEBInfo.end()) {
          ATH_MSG_DEBUG("Chain " << HLT::Identifier( chain ) << " adds PEBInfo " << it->second << " to stream " << st_type << "_" << st_name);
          // Copy ROB IDs directly
          robs = it->second.robs;
          // Convert SubDet IDs from uint32_t to eformat::SubDetector aka uint8_t
          for (const uint32_t subdetid : it->second.subdets) {
            subdets.insert(static_cast<eformat::SubDetector>( subdetid & 0xFF ));
          }
        }
      }

      eformat::helper::StreamTag streamTag(st_name, st_type, obeysLB, robs, subdets);
      ATH_CHECK(resultToFill.addStreamTag(streamTag));
    }
  }

  ATH_MSG_DEBUG("Number of streams for event " <<  resultToFill.getStreamTags().size() );
  return StatusCode::SUCCESS;
}

// =============================================================================

StatusCode StreamTagMakerTool::fillPEBInfoMap(std::unordered_map<DecisionID, PEBInfoWriterToolBase::PEBInfo>& map, const EventContext& ctx) const {
  for (const auto& key: m_pebDecisionKeys) {
    // Retrieve the decision container
    auto handle = SG::makeHandle(key, ctx);
    if (not handle.isValid())  {
      ATH_MSG_DEBUG("Missing input " <<  key.key() << " likely rejected");
      continue;
    }
    ATH_MSG_DEBUG("Processing decision container " << key.key());
    // Loop over decisions
    for (const Decision* d : *handle) {
      ATH_MSG_DEBUG("Processing decision " << *d);
      DecisionIDContainer ids;
      decisionIDs(d,ids);
      if (ids.empty()) {
        ATH_MSG_DEBUG("No chain passed for this decision object, skip retrieving PEB info");
        continue;
      }
      std::vector<uint32_t> robs;
      std::vector<uint32_t> subdets;
      if (d->getDetail(PEBInfoWriterToolBase::robListKey(), robs)) {
        ATH_MSG_DEBUG("Retrieved a list of " << robs.size() << " ROBs for this decision");
      }
      else {
        ATH_MSG_ERROR("Failed to retrieve " << PEBInfoWriterToolBase::robListKey() << " for decision container "
                      << key.key() << ", decision " << *d);
        return StatusCode::FAILURE;
      }
      if (d->getDetail(PEBInfoWriterToolBase::subDetListKey(), subdets)) {
        ATH_MSG_DEBUG("Retrieved a list of " << subdets.size() << " SubDets for this decision");
      }
      else {
        ATH_MSG_ERROR("Failed to retrieve " << PEBInfoWriterToolBase::subDetListKey() << " for decision container "
                      << key.key() << ", decision " << *d);
        return StatusCode::FAILURE;
      }
      // Assign PEBInfo to all passed chains for this decision
      for (const unsigned int id : ids) {
        // Convert leg id to chain id - changes nothing if id is already chain id
        const HLT::Identifier chainId = getIDFromLeg(id);
        ATH_MSG_DEBUG("Mapping PEBInfo to passed chain " << chainId.name());
        PEBInfoWriterToolBase::PEBInfo& pebi = map[chainId.numeric()];
        pebi.robs.insert(robs.begin(),robs.end());
        pebi.subdets.insert(subdets.begin(),subdets.end());
        if (pebi.robs.empty() && pebi.subdets.empty()) {
          // This would mean streaming the full event to a PEB stream
          ATH_MSG_ERROR("Empty PEB info for decision container " << key.key() << ", decision " << *d);
          return StatusCode::FAILURE;
        }
      }
    } // Loop over decisions
  } // Loop over decision containers
  return StatusCode::SUCCESS;
}
