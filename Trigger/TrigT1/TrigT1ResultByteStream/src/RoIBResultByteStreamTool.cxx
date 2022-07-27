/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// Local includes:
#include "RoIBResultByteStreamTool.h"

// Trigger includes:
#include "CTPfragment/CTPdataformat.h"
#include "CTPfragment/CTPfragment.h"
#include "L1TopoRDO/Helpers.h"
#include "L1TopoRDO/L1TopoRDO.h"
#include "TrigT1Result/L1TopoResult.h"
#include "TrigT1Result/RoIBResult.h"

// Athena includes
#include "CxxUtils/span.h"

// TDAQ includes:
#include "eformat/SourceIdentifier.h"

// System includes
#include <exception>
#include <sstream>

using DataType = OFFLINE_FRAGMENTS_NAMESPACE::PointerType;
using OFFLINE_FRAGMENTS_NAMESPACE::ROBFragment;

/**
 * The constructor takes care of correctly constructing the base class and
 * declaring the tool's interface to the framework.
 */
RoIBResultByteStreamTool::RoIBResultByteStreamTool( const std::string& type, const std::string& name,
                                                    const IInterface* parent )
  : base_class(type, name, parent) {}

/**
 * @brief Initialise the tool
 *
 * Fill the vector of configured ROB IDs and print the module IDs to debug message stream
 */
StatusCode RoIBResultByteStreamTool::initialize() {
  ATH_MSG_DEBUG("Initialising RoIBResultByteStreamTool");

  ConversionMode mode = getConversionMode(m_roibResultReadKey, m_roibResultWriteKey, msg());
  ATH_CHECK(mode!=ConversionMode::Undefined);

  // Set flags enabling/disabling each system
  m_doCTP = (m_ctpModuleID != 0xFF);
  m_doMuon = (m_muCTPIModuleID != 0xFF);
  m_doJetEnergy = !m_jetModuleID.empty();
  m_doEMTau = !m_emModuleID.empty();
  m_doTopo = !m_l1TopoModuleID.empty();

  if (msgLvl(MSG::DEBUG)) {
    std::string_view modeName {(mode==ConversionMode::Encoding) ? "encoding" : "decoding"};
    auto enabledStr = [](bool v) constexpr {return std::string_view(v ? "enabled" : "disabled");};
    ATH_MSG_DEBUG("CTP " << modeName << " is " << enabledStr(m_doCTP));
    ATH_MSG_DEBUG("MUCTPI " << modeName << " is " << enabledStr(m_doMuon));
    ATH_MSG_DEBUG("Jet/Energy " << modeName << " is " << enabledStr(m_doJetEnergy));
    ATH_MSG_DEBUG("EMTau " << modeName << " is " << enabledStr(m_doEMTau));
    ATH_MSG_DEBUG("L1Topo " << modeName << " is " << enabledStr(m_doTopo));
  }

  // Initialise data handles
  ATH_CHECK(m_roibResultWriteKey.initialize(mode==ConversionMode::Decoding));
  ATH_CHECK(m_roibResultReadKey.initialize(mode==ConversionMode::Encoding));
  ATH_CHECK(m_ctpRdoReadKey.initialize(m_doCTP && mode==ConversionMode::Encoding));

  // Set configured ROB IDs from module IDs
  std::vector<eformat::helper::SourceIdentifier> configuredROBSIDs;
  std::ostringstream str;
  ATH_MSG_DEBUG("Configured module IDs for:");

  // CTP
  if (m_doCTP) {
    configuredROBSIDs.emplace_back(eformat::TDAQ_CTP, m_ctpModuleID);
    ATH_MSG_DEBUG("   CTP                                  = 0x" << MSG::hex << m_ctpModuleID.value() << MSG::dec);
  }

  // MUCTPI
  if (m_doMuon) {
    configuredROBSIDs.emplace_back(eformat::TDAQ_MUON_CTP_INTERFACE, m_muCTPIModuleID);
    ATH_MSG_DEBUG("   muCTPi                               = 0x" << MSG::hex << m_muCTPIModuleID.value() << MSG::dec);
  }

  // Jet/Energy
  if (m_doJetEnergy) {
    str.str("");
    for (const uint16_t module_id : m_jetModuleID) {
      configuredROBSIDs.emplace_back(eformat::TDAQ_CALO_JET_PROC_ROI, module_id);
      str << "0x" << std::hex << module_id << std::dec << " ";
    }
    ATH_MSG_DEBUG("   Calorimeter Jet/Energy Processor RoI = " << str.str());
  }

  // EM/Tau
  if (m_doEMTau) {
    str.str("");
    for (const uint16_t module_id : m_emModuleID) {
      configuredROBSIDs.emplace_back(eformat::TDAQ_CALO_CLUSTER_PROC_ROI, module_id);
      str << "0x" << std::hex << module_id << std::dec << " ";
    }
    ATH_MSG_DEBUG("   Calorimeter Cluster Processor RoI    = " << str.str());
  }

  // L1Topo
  if (m_doTopo) {
    str.str("");
    for (const uint16_t module_id : m_l1TopoModuleID) {
      configuredROBSIDs.emplace_back(eformat::TDAQ_CALO_TOPO_PROC, module_id);
      str << "0x" << std::hex << module_id << std::dec << " ";
    }
    ATH_MSG_DEBUG("   L1Topo                               = " << str.str());
  }

  // Fill the ROB ID vector
  for (const auto& sid : configuredROBSIDs) {
    m_configuredROBIds.push_back( sid.code() );
  }

  // Retrieve the ByteStreamCnvSvc for updating event header when writing CTP data to BS
  if (m_doCTP && mode==ConversionMode::Encoding) {
    m_byteStreamEventAccess = serviceLocator()->service("ByteStreamCnvSvc");
    ATH_CHECK(m_byteStreamEventAccess.isValid());
  }

  return StatusCode::SUCCESS;
}

StatusCode RoIBResultByteStreamTool::convertToBS(std::vector<OFFLINE_FRAGMENTS_NAMESPACE_WRITE::ROBFragment*>& vrobf,
                                                 const EventContext& eventContext) {
  auto roibResult = SG::makeHandle(m_roibResultReadKey, eventContext);
  ATH_CHECK(roibResult.isValid());
  ATH_MSG_DEBUG("Obtained ROIB::RoIBResult with key " << m_roibResultReadKey.key() << " for conversion to ByteStream");

  auto addRob = [&](const eformat::helper::SourceIdentifier& sid, const size_t ndata, const uint32_t* data){
    vrobf.push_back(newRobFragment(eventContext, sid.code(), ndata, data, m_detEvType));
    return vrobf.back();
  };
  auto convertDataToRob = [&](const eformat::helper::SourceIdentifier& sid, const auto& dataVec){
    uint32_t* data = newRodData(eventContext, dataVec.size());
    ATH_MSG_VERBOSE("Dumping words for " << sid.human());
    for (size_t i=0; i<dataVec.size(); ++i) {
      ATH_MSG_VERBOSE("  0x" << MSG::hex << std::setw(8) << dataVec.at(i).roIWord() << MSG::dec);
      data[i] = dataVec.at(i).roIWord();
    }
    return addRob(sid, dataVec.size(), data);
  };

  // CTP
  if (m_doCTP) {
    const eformat::helper::SourceIdentifier ctpSID{eformat::TDAQ_CTP, m_ctpModuleID};
    const std::vector<ROIB::CTPRoI>& ctpDataVec = roibResult->cTPResult().roIVec();
    OFFLINE_FRAGMENTS_NAMESPACE_WRITE::ROBFragment* rob = convertDataToRob(ctpSID, ctpDataVec);

    // Update ROD minor version in the ROB based on information from CTP_RDO
    auto ctpRdo = SG::makeHandle(m_ctpRdoReadKey, eventContext);
    ATH_CHECK(ctpRdo.isValid());
    const uint16_t rodMinorVersion = ctpRodMinorVersion(*ctpRdo);
    rob->rod_minor_version(rodMinorVersion);

    // Update L1 bits in event header
    // TODO: change to context-aware method when it is added to the IByteStreamEventAccess interface
    RawEventWrite* rawEvent = m_byteStreamEventAccess->getRawEvent();
    std::bitset<512> tbp = ROIB::convertToBitset(roibResult->cTPResult().TBP());
    std::bitset<512> tap = ROIB::convertToBitset(roibResult->cTPResult().TAP());
    std::bitset<512> tav = ROIB::convertToBitset(roibResult->cTPResult().TAV());
    constexpr size_t word_size = 32; // ATLAS raw event format word size
    constexpr size_t words_per_set = 512 / word_size; // 512 CTP bits / word size
    constexpr size_t num_words = 3 * words_per_set; // 3 sets: TBP, TAP, TAV
    uint32_t* l1bits_data = newRodData(eventContext, num_words);
    size_t iset{0};
    for (const std::bitset<512>& bset : {tbp, tap, tav}) {
      const std::string sbits = bset.to_string();
      const size_t iword_output_start = words_per_set * iset;
      for (size_t iword_in_set = 0; iword_in_set < words_per_set; ++iword_in_set) {
        const size_t bword_pos = 512 - (iword_in_set+1) * word_size;
        std::bitset<word_size> bword{sbits.substr(bword_pos, word_size)};
        l1bits_data[iword_output_start + iword_in_set] = static_cast<uint32_t>(bword.to_ulong());
      }
      ++iset;
    }
    rawEvent->lvl1_trigger_info(num_words, l1bits_data);

    // Update L1 TriggerType in event header
    const uint32_t triggerType = roibResult->cTPResult().header().triggerType();
    rawEvent->lvl1_trigger_type(static_cast<uint8_t>(triggerType & 0xFF));
  }

  // Muon
  if (m_doMuon) {
    const eformat::helper::SourceIdentifier muctpiSID{eformat::TDAQ_MUON_CTP_INTERFACE, m_muCTPIModuleID};
    const std::vector<ROIB::MuCTPIRoI>& muctpiDataVec = roibResult->muCTPIResult().roIVec();
    convertDataToRob(muctpiSID, muctpiDataVec);
  }

  // Jet/Energy
  if (m_doJetEnergy) {
    const std::vector<ROIB::JetEnergyResult>& jetEnergyResultVec = roibResult->jetEnergyResult();
    for (size_t slink=0; slink<jetEnergyResultVec.size(); ++slink) {
      const eformat::helper::SourceIdentifier jetSID{eformat::TDAQ_CALO_JET_PROC_ROI, m_jetModuleID.value().at(slink)};
      convertDataToRob(jetSID, jetEnergyResultVec.at(slink).roIVec());
    }
  }

  // EMTau
  if (m_doEMTau) {
    const std::vector<ROIB::EMTauResult>& emTauResultVec = roibResult->eMTauResult();
    for (size_t slink=0; slink<emTauResultVec.size(); ++slink) {
      const eformat::helper::SourceIdentifier emSID{eformat::TDAQ_CALO_CLUSTER_PROC_ROI, m_emModuleID.value().at(slink)};
      convertDataToRob(emSID, emTauResultVec.at(slink).roIVec());
    }
  }

  // L1Topo (slightly different interface from the others)
  if (m_doTopo) {
    const std::vector<ROIB::L1TopoResult>& l1TopoResultVec = roibResult->l1TopoResult();
    for (size_t slink=0; slink<l1TopoResultVec.size(); ++slink) {
      eformat::helper::SourceIdentifier topoSID{l1TopoResultVec.at(slink).rdo().getSourceID()};
      ATH_MSG_VERBOSE("L1Topo source ID from RDO "  << L1Topo::formatHex8(topoSID.code()));
      if (topoSID.code() == 0 && slink < m_l1TopoModuleID.size()) {
        topoSID = eformat::helper::SourceIdentifier(eformat::TDAQ_CALO_TOPO_PROC, m_l1TopoModuleID.value().at(slink));
        ATH_MSG_DEBUG("Source ID in L1TopoRDO was zero so using Property for slink " << slink << ": "
                      << L1Topo::formatHex8(topoSID.code()));
      }
      else if (topoSID.code() == 0) {
        topoSID = eformat::helper::SourceIdentifier( eformat::TDAQ_CALO_TOPO_PROC, 0 );
        ATH_MSG_WARNING("Source ID in L1TopoRDO was zero, no properties available for slink counter " << slink
                        << ", so as a fall back, constructed module 0 with source ID "
                        << L1Topo::formatHex8(topoSID.code()));
      }
      const std::vector<uint32_t>& dataVec = l1TopoResultVec.at(slink).rdo().getDataWords();
      uint32_t* data = newRodData(eventContext, dataVec.size());
      for (size_t i=0; i<dataVec.size(); ++i) {
        ATH_MSG_VERBOSE("  " << MSG::hex << std::setw(8) << std::showbase << dataVec.at(i) << std::noshowbase << MSG::dec);
        data[i] = dataVec.at(i);
      }
      addRob(topoSID, dataVec.size(), data);
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode RoIBResultByteStreamTool::convertFromBS(const std::vector<const ROBFragment*>& vrobf,
                                                   const EventContext& eventContext) const {

  if (m_roibResultWriteKey.empty()) {
    ATH_MSG_ERROR("Conversion from BS to xAOD RoI requested but RoI WriteHandleKey is empty");
    return StatusCode::FAILURE;
  }
  // Create all RDOs
  ROIB::CTPResult cTPResult;
  ROIB::MuCTPIResult muCTPIResult;
  std::vector< ROIB::JetEnergyResult > jetEnergyResult(2);
  std::vector< ROIB::EMTauResult > eMTauResult(4);
  std::vector< ROIB::L1TopoResult > l1TopoResult;

  // Create flags indicating whether or not ROB fragment was found
  bool cTPFound = false;
  bool muCTPIFound = false;
  bool jetEnergyFound[2] = {false, false};
  bool eMTauFound[4] = {false, false, false, false};
  bool l1TopoFound = false;

  // Loop over ROB fragments
  for (const ROBFragment* p_robf : vrobf) {
    const ROBFragment& robf = *p_robf;
    eformat::helper::SourceIdentifier rodSID(robf.rod_source_id());
    eformat::SubDetector rodSubdetId = rodSID.subdetector_id();
    uint16_t rodModuleId = rodSID.module_id();

    switch (rodSubdetId) {
      // -----------------------------------------------------------------------
      // CTP
      // -----------------------------------------------------------------------
      case eformat::TDAQ_CTP: {
        if (!m_doCTP) continue;

        // Check if the current ROD module ID matches the configured CTP module ID.
        // Accept also ctpModuleId=0 to catch the early data with the old CTP firmware, which assigns 0x770000 to both
        // DAQ and LVL2 ROD fragment - for this data we readout always multiple bunches for the DAQ ROD, therefore
        // rodFragSize=46 should identify the LVL2 ROD fragment
        if (rodModuleId != m_ctpModuleID && rodModuleId != 0) continue;
        if (rodModuleId == 0 && robf.rod_fragment_size_word() != 46) continue;

        // Flag as found
        cTPFound = true;
        ATH_MSG_DEBUG("   Found CTP ROD with source ID " << MSG::hex << rodSID.code() << MSG::dec);

        // Extract the data
        DataStatus status;
        ROIB::Header header = roibHeader(robf, status);
        std::vector<ROIB::CTPRoI> content = roibContent<ROIB::CTPRoI>(robf);
        ROIB::Trailer trailer = roibTrailer(status, content.size());

        // Extract the CTP version number
        const uint32_t* rod;
        robf.rod_start(rod);
        unsigned int ctpVersionNumber = ((rod[CTPdataformat::Helper::FormatVersionPos] >> CTPdataformat::CTPFormatVersionShift) & CTPdataformat::CTPFormatVersionMask);

        // Create CTPResult object
        cTPResult = ROIB::CTPResult(ctpVersionNumber, std::move(header), std::move(trailer), std::move(content));
        break;
      }

      // -----------------------------------------------------------------------
      // MUCTPI
      // -----------------------------------------------------------------------
      case eformat::TDAQ_MUON_CTP_INTERFACE: {
        if (!m_doMuon) continue;

        // Check if the current ROD module ID matches the configured MuCTPI module ID
        if (rodModuleId != m_muCTPIModuleID) continue;

        // Flag as found
        muCTPIFound = true;
        ATH_MSG_DEBUG("   Found MuCTPI ROD with source ID " << MSG::hex << rodSID.code() << MSG::dec);

        // Extract the data
        DataStatus status;
        ROIB::Header header = roibHeader(robf, status);
        std::vector<ROIB::MuCTPIRoI> content = roibContent<ROIB::MuCTPIRoI>(robf);
        ROIB::Trailer trailer = roibTrailer(status, content.size());

        // Create MuCTPIResult object
        muCTPIResult = ROIB::MuCTPIResult(std::move(header), std::move(trailer), std::move(content));
        break;
      }

      // -----------------------------------------------------------------------
      // Jet/Energy
      // -----------------------------------------------------------------------
      case eformat::TDAQ_CALO_JET_PROC_ROI: {
        if (!m_doJetEnergy) continue;

        // Check if the current ROD module ID matches the configured Jet/Energy module IDs
        auto it = std::find(m_jetModuleID.begin(), m_jetModuleID.end(), rodModuleId);
        if (it == m_jetModuleID.end()) continue;
        size_t index = static_cast<size_t>(std::distance(m_jetModuleID.begin(), it));

        // Flag as found
        jetEnergyFound[index] = true;
        ATH_MSG_DEBUG("   Found Jet/Energy ROD with source ID " << MSG::hex << rodSID.code() << MSG::dec);

        // Extract the data
        DataStatus status;
        ROIB::Header header = roibHeader(robf, status);
        std::vector<ROIB::JetEnergyRoI> content = roibContent<ROIB::JetEnergyRoI>(robf);
        ROIB::Trailer trailer = roibTrailer(status, content.size());

        // Create JetEnergyResult object
        jetEnergyResult[index] = ROIB::JetEnergyResult(std::move(header), std::move(trailer), std::move(content));
        break;
      }

      // -----------------------------------------------------------------------
      // EM/Tau
      // -----------------------------------------------------------------------
      case eformat::TDAQ_CALO_CLUSTER_PROC_ROI: {
        if (!m_doEMTau) continue;

        // Check if the current ROD module ID matches the configured EM/Tau module IDs
        auto it = std::find(m_emModuleID.begin(), m_emModuleID.end(), rodModuleId);
        if (it == m_emModuleID.end()) continue;
        size_t index = static_cast<size_t>(std::distance(m_emModuleID.begin(), it));

        // Flag as found
        eMTauFound[index] = true;
        ATH_MSG_DEBUG("   Found EM/Tau ROD with source ID " << MSG::hex << rodSID.code() << MSG::dec);

        // Extract the data
        DataStatus status;
        ROIB::Header header = roibHeader(robf, status);
        std::vector<ROIB::EMTauRoI> content = roibContent<ROIB::EMTauRoI>(robf);
        ROIB::Trailer trailer = roibTrailer(status, content.size());

        // Create EMTauResult object
        eMTauResult[index] = ROIB::EMTauResult(std::move(header), std::move(trailer), std::move(content));
        break;
      }

      // -----------------------------------------------------------------------
      // L1Topo
      // -----------------------------------------------------------------------
      case eformat::TDAQ_CALO_TOPO_PROC: {
        if (!m_doTopo) continue;

        // Check if the current ROD module ID matches the configured L1Topo module IDs
        auto it = std::find(m_l1TopoModuleID.begin(), m_l1TopoModuleID.end(), rodModuleId);
        if (it == m_l1TopoModuleID.end()) continue;

        // Flag as found
        l1TopoFound = true;
        ATH_MSG_DEBUG("   Found L1Topo ROD with source ID " << MSG::hex << rodSID.code() << MSG::dec);

        // Extract the data
        DataStatus status;
        ROIB::Header header = roibHeader(robf, status);
        L1TopoRDO content = l1topoContent(robf);
        ROIB::Trailer trailer = roibTrailer(status, content.getDataWords().size());

        // Set status words
        content.setStatusWords({status.status_word, status.status_info});

        // Flag errors in RDO
        if (status.status_word != 0) content.setError(L1Topo::Error::SLINK_STATUS_ERROR);
        if (status.rob_error) content.setError(L1Topo::Error::ROB_ERROR);
        if (status.rod_error) content.setError(L1Topo::Error::ROD_ERROR);

        // Create L1TopoResult object
        l1TopoResult.emplace_back(std::move(header), std::move(trailer), std::move(content));
        break;
      }

      default: {
        ATH_MSG_DEBUG("Skipping ROD with SubDetID " << rodSID.human_detector());
        break;
      }
    }

  } // End of loop over all ROB fragments

  ATH_MSG_DEBUG("Building RoIBResult with the following inputs:");
  ATH_MSG_DEBUG("  CTP             - " << (cTPFound ? "found" : "not found"));
  ATH_MSG_DEBUG("  MUCTPI          - " << (muCTPIFound ? "found" : "not found"));
  ATH_MSG_DEBUG("  Jet/Energy[0/1] - " << (jetEnergyFound[0] ? "found" : "not found") << "/"
                                       << (jetEnergyFound[1] ? "found" : "not found"));
  ATH_MSG_DEBUG("  EM/Tau[0/1/2/3] - " << (eMTauFound[0] ? "found" : "not found") << "/"
                                       << (eMTauFound[1] ? "found" : "not found") << "/"
                                       << (eMTauFound[2] ? "found" : "not found") << "/"
                                       << (eMTauFound[3] ? "found" : "not found"));
  ATH_MSG_DEBUG("  L1Topo          - " << (l1TopoFound ? "found" : "not found"));


  // Create and record the RoIBResult
  auto roibResult = SG::makeHandle(m_roibResultWriteKey, eventContext);
  ATH_CHECK(roibResult.record(std::make_unique<ROIB::RoIBResult>(std::move(muCTPIResult), std::move(cTPResult), std::move(jetEnergyResult), std::move(eMTauResult))));
  ATH_MSG_DEBUG("Recorded RoIBResult with key " << m_roibResultWriteKey.key());

  // Add L1Topo result
  if (l1TopoFound) roibResult->l1TopoResult(std::move(l1TopoResult));

  return StatusCode::SUCCESS;
}

ROIB::Header RoIBResultByteStreamTool::roibHeader(const ROBFragment& rob,
                                                  DataStatus& dataStatus) const {
  // Get format version and event number of fragment
  uint32_t formatVersion = rob.rod_version();
  uint32_t evtNum = rob.rod_lvl1_id();
  uint32_t robFragSize = rob.fragment_size_word();
  uint32_t rodFragSize = rob.rod_fragment_size_word();
  uint32_t robId = rob.source_id();
  uint32_t rodId = rob.rod_source_id();

  // Check for errors
  dataStatus.rob_error = false;
  dataStatus.rod_error = false;
  dataStatus.status_word = 0;
  dataStatus.status_info = 0;

  try {
    if (rob.check_rob()) ATH_MSG_VERBOSE("ROB fragment checked ok");
  }
  catch (std::exception const & ex) {
    ATH_MSG_WARNING("ROB fragment not valid: " << ex.what());
    dataStatus.rob_error = true;
  }

  try {
    if (rob.check_rod()) ATH_MSG_VERBOSE("ROD fragment checked ok");
  }
  catch (std::exception const & ex) {
    ATH_MSG_WARNING("ROD fragment not valid: " << ex.what());
    dataStatus.rod_error = true;
  }

  // Process the status words
  DataType status;
  rob.rod_status(status);
  uint32_t nstatus = rob.rod_nstatus();
  ATH_MSG_VERBOSE("Number of status words: " << nstatus);
  if (msgLvl(MSG::VERBOSE)) {
    for (uint32_t i=0; i<nstatus; ++i, ++status) {
      ATH_MSG_VERBOSE("   Status word: 0x" << MSG::hex << std::setw(8) << std::setfill('0') << *status);
    }
  }
  rob.rod_status(status);

  if(nstatus > 0) dataStatus.status_word = static_cast<uint32_t>(*status);
  if(nstatus > 1) {
      ++status;
      dataStatus.status_info = static_cast<uint32_t>(*status);
  }

  ATH_MSG_DEBUG("ROB ID 0x" << MSG::hex << robId <<  " ROD ID 0x" << rodId << MSG::dec << " ROB fragment size "
                << robFragSize << " ROD fragment size " << rodFragSize);

  return ROIB::Header(rodId, evtNum, formatVersion);
}

ROIB::Trailer RoIBResultByteStreamTool::roibTrailer(const DataStatus& dataStatus, uint32_t dataSize) const {
  std::vector<uint32_t> words;
  words.reserve(5);
  words.push_back(dataStatus.status_word); // error status
  words.push_back(dataStatus.status_info); // status info
  words.push_back(2);                      // number of status words
  words.push_back(dataSize);               // number of data words
  words.push_back(1);                      // status block position
  return ROIB::Trailer(std::move(words));
}

template<typename RoIType> std::vector<RoIType> RoIBResultByteStreamTool::roibContent(const ROBFragment& rob) const {
  uint32_t ndata = rob.rod_ndata();
  DataType data;
  rob.rod_data(data);
  std::vector<RoIType> content;
  content.reserve(ndata);
  ATH_MSG_VERBOSE("   Dumping RoI Words:");
  for (const uint32_t word : CxxUtils::span{data, ndata}) {
    ATH_MSG_VERBOSE("       0x" << MSG::hex << std::setfill('0') << std::setw(8) << word << MSG::dec);
    content.push_back(RoIType{word});
  }
  return content;
}

L1TopoRDO RoIBResultByteStreamTool::l1topoContent(const ROBFragment& rob) const {
  uint32_t ndata = rob.rod_ndata();
  DataType data;
  rob.rod_data(data);
  L1TopoRDO content;
  ATH_MSG_VERBOSE( "   Dumping RoI Words:" );
  std::vector<uint32_t> vDataWords;
  vDataWords.reserve(ndata);
  for (const uint32_t word : CxxUtils::span{data, ndata}) {
    ATH_MSG_VERBOSE("       0x" << MSG::hex << std::setfill('0') << std::setw(8) << word << MSG::dec);
    vDataWords.push_back(word);
  }
  content.setDataWords(std::move(vDataWords));
  content.setSourceID(rob.rod_source_id());
  return content;
}

uint16_t RoIBResultByteStreamTool::ctpRodMinorVersion(const CTP_RDO& rdo) const {
  const CTPdataformatVersion& ctpFormat = rdo.getCTPVersion();
  const uint32_t ctpVersionNumber = rdo.getCTPVersionNumber();
  const uint32_t l1aPosition = rdo.getL1AcceptBunchPosition();
  const uint32_t nAddWords = rdo.getNumberOfAdditionalWords();

  ATH_MSG_DEBUG("CTP version number: " << ctpVersionNumber);
  ATH_MSG_DEBUG("L1A bunch position: " << l1aPosition);
  ATH_MSG_DEBUG("N additional words: " << nAddWords);
  ATH_MSG_DEBUG("N bunches: " << rdo.getNumberOfBunches());

  uint16_t rodMinorVersion{0};
  rodMinorVersion |= ((nAddWords & ctpFormat.getProgrammableExtraWordsMask()) << ctpFormat.getProgrammableExtraWordsShift());
  rodMinorVersion |= ((l1aPosition & ctpFormat.getL1APositionMask()) << ctpFormat.getL1APositionShift());
  rodMinorVersion |= ((ctpVersionNumber & ctpFormat.getCTPFormatVersionMask()) << ctpFormat.getCTPFormatVersionShift());

  ATH_MSG_DEBUG("CTP ROD minor version calculated as 0x" << MSG::hex << rodMinorVersion << MSG::dec);

  return rodMinorVersion;
}
