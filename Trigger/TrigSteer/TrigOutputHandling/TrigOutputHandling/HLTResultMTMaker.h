/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGOUTPUTHANDLING_HLTRESULTMTMAKER_H
#define TRIGOUTPUTHANDLING_HLTRESULTMTMAKER_H

// Trigger includes
#include "TrigSteeringEvent/HLTResultMT.h"
#include "TrigOutputHandling/HLTResultMTMakerTool.h"

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaMonitoringKernel/GenericMonitoringTool.h"
#include "ByteStreamData/ByteStreamMetadataContainer.h"
#include "StoreGate/ReadHandle.h"
#include "StoreGate/WriteHandle.h"
#include "xAODTrigger/TrigCompositeContainer.h"

// Gaudi includes
#include "Gaudi/Interfaces/IOptionsSvc.h"

// System includes
#include <atomic>

/** @class HLTResultMTMaker
 *  @brief Tool to create the HLTResultMT at the end of each event
 **/
class HLTResultMTMaker : public AthAlgTool {
public:
  /// Standard constructor
  HLTResultMTMaker(const std::string& type, const std::string& name, const IInterface* parent);
  /// Standard destructor
  virtual ~HLTResultMTMaker() override = default;

  // ------------------------- IStateful methods -------------------------------
  virtual StatusCode initialize() override;
  virtual StatusCode start() override;
  virtual StatusCode finalize() override;

  // ------------------------- Specific methods of this tool -------------------
  /// Create and fill a new HLTResultMT, and record it in the event store
  StatusCode makeResult(const EventContext& eventContext) const;
  /// Fill an existing HLTResultMT with event information
  StatusCode fillResult(HLT::HLTResultMT& hltResult, const EventContext& eventContext) const;
  /// Return name of the HLTResultMT
  const std::string& resultName() const {return m_hltResultWHKey.key();}

private:
  // ------------------------- Methods -----------------------------------------
  /// Check ROB and SubDet lists in StreamTags and remove those which are disabled
  void validatePEBInfo(HLT::HLTResultMT& hltResult) const;
  /// Check ExtraROBs and ExtraSubDets against the detector mask and drop the masked out IDs
  StatusCode validateExtraROBsAndSubDets();

  // ------------------------- Properties --------------------------------------
  /// StoreGate key for the HLTResultMT
  SG::WriteHandleKey<HLT::HLTResultMT> m_hltResultWHKey {
    this, "HLTResultWHKey", "HLTResultMT",
    "Key of the output HLTResultMT object"
  };
  /// StoreGate key for HLT Runtime Metadata container
  SG::WriteHandleKey<xAOD::TrigCompositeContainer> m_runtimeMetadataWHKey { 
    this, "HLTRuntimeMetadataWHKey", "HLT_RuntimeMetadata",
    "Key of the output with additional runtime metadata" };
  /// StoreGate key for the ByteStreamMetadata container to retrieve detector mask
  SG::ReadHandleKey<ByteStreamMetadataContainer> m_bsMetaDataContRHKey {
    this, "ByteStreamMetadataRHKey", "InputMetaDataStore+ByteStreamMetadata",
    "Key of the ByteStreamMetadataContainer to retrieve the detector mask"
  };
  /// Tool creating stream tags (defines if event is accepted)
  ToolHandle<HLTResultMTMakerTool> m_streamTagMaker {
    this, "StreamTagMaker", "",
    "Tool creating stream tags (defines if event is accepted)"
  };
  /// Tools filling the HLTResultMT object
  ToolHandleArray<HLTResultMTMakerTool> m_makerTools {
    this, "MakerTools", {},
    "Set of additional tools that fill content of the HLTResultMT"
  };
  /// Monitoring tool
  ToolHandle<GenericMonitoringTool> m_monTool {
    this, "MonTool", "",
    "Monitoring tool"
  };
  /// Handle to JobOptionsSvc used to retrieve the DataFlowConfig property
  ServiceHandle<Gaudi::Interfaces::IOptionsSvc> m_jobOptionsSvc {
    this, "JobOptionsSvc", "JobOptionsSvc",
    "Job options service to retrieve DataFlowConfig"
  };
  /// Extra enabled ROBs
  Gaudi::Property<std::vector<uint32_t>> m_extraEnabledROBs {
    this, "ExtraEnabledROBs", {},
    "Extra ROBs which can be requested in a stream tag, but are not part of the ROS-ROB map"
  };
  /// Extra enabled ROBs checked against detector mask
  Gaudi::Property<std::vector<uint32_t>> m_extraROBs {
    this, "ExtraROBs", {},
    "Same as ExtraEnabledROBs but checked against detector mask and dropped if SubDet is masked out"
  };
  /// Extra enabled SubDets
  Gaudi::Property<std::vector<uint32_t>> m_extraEnabledSubDets {
    this, "ExtraEnabledSubDets", {},
    "Extra SubDets which can be requested in a stream tag, but are not part of the ROS-ROB map"
  };
  /// Extra enabled SubDets checked against detector mask
  Gaudi::Property<std::vector<uint32_t>> m_extraSubDets {
    this, "ExtraSubDets", {},
    "Same as ExtraEnabledSubDets but checked against detector mask and dropped if SubDet is masked out"
  };

  // ------------------------- Other private members ---------------------------
  /// List of enabled ROBs retrieved during initialisation
  std::set<uint32_t> m_enabledROBs;
  /// List of enabled SubDets retrieved during initialisation
  std::set<eformat::SubDetector> m_enabledSubDets;
  /// If true, don't call validatePEBInfo
  bool m_skipValidatePEBInfo {false};
  /// Flag if empty PEB list error was already printed
  mutable std::atomic_bool m_emptyPEBInfoErrorPrinted{false};
};

#endif // TRIGOUTPUTHANDLING_HLTRESULTMTMAKER_H
