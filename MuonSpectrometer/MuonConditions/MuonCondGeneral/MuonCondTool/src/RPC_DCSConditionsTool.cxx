/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "RPC_DCSConditionsTool.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>

#include "AthenaPoolUtilities/AthenaAttributeList.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "CoralBase/Attribute.h"
#include "CoralBase/AttributeListSpecification.h"
#include "GaudiKernel/MsgStream.h"
#include "Identifier/Identifier.h"
#include "Identifier/IdentifierHash.h"
#include "MuonCondSvc/MdtStringUtils.h"
#include "MuonIdHelpers/RpcIdHelper.h"
#include "PathResolver/PathResolver.h"
#include "SGTools/TransientAddress.h"

//**********************************************************
//* Author Monica Verducci monica.verducci@cern.ch
//*
//* Tool to retrieve the RPC DCS Info from COOL DB
//* retrieving of tables from DB
//*********************************************************

RPC_DCSConditionsTool::RPC_DCSConditionsTool(const std::string& type, const std::string& name, const IInterface* parent) :
    AthAlgTool(type, name, parent),
    m_IOVSvc(nullptr),
    m_rpcIdHelper(nullptr),
    m_log(msgSvc(), name),
    m_debug(false),
    m_verbose(false),
    m_DataLocation("keyRPCDCS"),
    m_chronoSvc(nullptr) {
    declareInterface<IRPC_DCSConditionsTool>(this);

    declareProperty("OffPanelFolder", m_offPanelFolder = "/RPC/DCS/OffRopanels");
    declareProperty("DeadPanel", m_deadPanelFolder = "/RPC/DCS/DeadRopanels");
    m_RPCPaneloff.str("EMPTY");
    m_RPCPaneldead.str("EMPTY");
}

// StatusCode RPC_DCSConditionsTool::updateAddress(SG::TransientAddress* /*tad*/)
StatusCode RPC_DCSConditionsTool::updateAddress(StoreID::type /*storeID*/, SG::TransientAddress* /*tad*/, const EventContext& /*ctx*/) {
    return StatusCode::FAILURE;
}

StatusCode RPC_DCSConditionsTool::initialize() {
    m_log.setLevel(msgLevel());
    m_debug = m_log.level() <= MSG::DEBUG;
    m_verbose = m_log.level() <= MSG::VERBOSE;

    m_log << MSG::INFO << "Initializing - folders names are: Panel Off  " << m_offPanelFolder << " / Panel Dead" << m_deadPanelFolder
          << endmsg;

    StatusCode sc = detStore()->retrieve(m_rpcIdHelper, "RPCIDHELPER");
    if (sc.isFailure()) {
        m_log << MSG::FATAL << " Cannot retrieve RpcIdHelper " << endmsg;
        return sc;
    }

    // Get interface to IOVSvc
    m_IOVSvc = nullptr;
    bool CREATEIF(true);
    sc = service("IOVSvc", m_IOVSvc, CREATEIF);
    if (sc.isFailure()) {
        m_log << MSG::ERROR << "Unable to get the IOVSvc" << endmsg;
        return StatusCode::FAILURE;
    }

    if (sc.isFailure()) return StatusCode::FAILURE;

    // initialize the chrono service
    sc = service("ChronoStatSvc", m_chronoSvc);
    if (sc != StatusCode::SUCCESS) {
        m_log << MSG::ERROR << "Could not find the ChronoSvc" << endmsg;
        return sc;
    }

    if (sc.isFailure()) return StatusCode::FAILURE;

    return StatusCode::SUCCESS;
}

StatusCode RPC_DCSConditionsTool::loadParameters(IOVSVC_CALLBACK_ARGS_P(I, keys)) {
    m_log.setLevel(msgLevel());
    m_debug = m_log.level() <= MSG::DEBUG;
    m_verbose = m_log.level() <= MSG::VERBOSE;

    std::list<std::string>::const_iterator itr;
    for (itr = keys.begin(); itr != keys.end(); ++itr) {
        m_log << MSG::INFO << "LoadParameters Dead and Off" << *itr << " I=" << I << " " << endmsg;
        if (*itr == m_deadPanelFolder) {
            StatusCode sc = loadPanelDead(I, keys);
            if (sc.isFailure()) { return sc; }
        } else if (*itr == m_offPanelFolder) {
            StatusCode sc = loadPanelOff(I, keys);
            if (sc.isFailure()) { return sc; }
        }
    }

    return StatusCode::SUCCESS;
}

StatusCode RPC_DCSConditionsTool::loadPanelOff(IOVSVC_CALLBACK_ARGS_P(I, keys)) {
    m_log.setLevel(msgLevel());
    m_debug = m_log.level() <= MSG::DEBUG;
    m_verbose = m_log.level() <= MSG::VERBOSE;

    StatusCode sc = StatusCode::SUCCESS;
    m_log << MSG::INFO << "Load Off Panel from DCS DB" << endmsg;

    // Print out callback information
    if (m_debug) m_log << MSG::DEBUG << "Level " << I << " Keys: ";
    std::list<std::string>::const_iterator keyIt = keys.begin();
    for (; keyIt != keys.end(); ++keyIt) m_log << MSG::DEBUG << *keyIt << " ";
    if (m_debug) m_log << MSG::DEBUG << endmsg;

    const CondAttrListCollection* atrc;
    m_log << MSG::INFO << "Try to read from folder <" << m_offPanelFolder << ">" << endmsg;

    sc = detStore()->retrieve(atrc, m_offPanelFolder);
    if (sc.isFailure()) {
        m_log << MSG::ERROR << "could not retreive the CondAttrListCollection from DB folder " << m_offPanelFolder << endmsg;
        return sc;
    }

    else if (m_debug)
        m_log << MSG::DEBUG << " CondAttrListCollection from DB folder have been obtained with size " << atrc->size() << endmsg;

    CondAttrListCollection::const_iterator itr;
    for (itr = atrc->begin(); itr != atrc->end(); ++itr) {
        // itr=atrc->chanAttrListPair(chanNum);
        const coral::AttributeList& atr = itr->second;

        std::string panel_off;
        std::string panel_reason_off;
        // if(atr.size()==1){
        if (atr.size()) {
            panel_off = *(static_cast<const std::string*>((atr["RpcOffROPanelIds"]).addressOfData()));
            panel_reason_off = *(static_cast<const std::string*>((atr["RpcOffROPanelReasons"]).addressOfData()));

            if (m_debug) m_log << MSG::DEBUG << "panel_off " << panel_off << endmsg;
            if (m_debug) m_log << MSG::DEBUG << "panel_reason " << panel_reason_off << endmsg;

            char delimiter = ',';
            const auto info_panel = MuonCalib::MdtStringUtils::tokenize(panel_off, delimiter);

            Identifier PanelId;

            for (unsigned int i = 0; i < info_panel.size(); i++) {
                const std::string_view& ch_tmp = info_panel[i];
                if (m_debug) m_log << MSG::DEBUG << " info_panel " << ch_tmp << " " << MuonCalib::MdtStringUtils::atoi(ch_tmp) << endmsg;

                PanelId = MuonCalib::MdtStringUtils::atoi(ch_tmp);

                if (PanelId.get_compact()) {
                    if (m_debug) m_log << MSG::DEBUG << "OFFPANEL " << m_rpcIdHelper->show_to_string(PanelId) << endmsg;
                    Identifier atlasId = m_rpcIdHelper->panelID(PanelId);
                    if (atlasId != 0) m_cachedOffPanelId.push_back(PanelId);
                    if (m_debug) m_log << MSG::DEBUG << "push-back" << endmsg;
                }
            }
        }
    }

    return StatusCode::SUCCESS;
}

StatusCode RPC_DCSConditionsTool::loadPanelDead(IOVSVC_CALLBACK_ARGS_P(I, keys)) {
    m_log.setLevel(msgLevel());
    m_debug = m_log.level() <= MSG::DEBUG;
    m_verbose = m_log.level() <= MSG::VERBOSE;

    StatusCode sc = StatusCode::SUCCESS;
    m_log << MSG::INFO << "Load Dead Panel from DCS DB" << endmsg;

    // Print out callback information
    if (m_debug) m_log << MSG::DEBUG << "Level " << I << " Keys: ";
    std::list<std::string>::const_iterator keyIt = keys.begin();
    for (; keyIt != keys.end(); ++keyIt)
        if (m_debug) m_log << MSG::DEBUG << *keyIt << " ";
    if (m_debug) m_log << MSG::DEBUG << endmsg;

    const CondAttrListCollection* atrc;
    m_log << MSG::INFO << "Try to read from folder <" << m_deadPanelFolder << ">" << endmsg;

    sc = detStore()->retrieve(atrc, m_deadPanelFolder);
    if (sc.isFailure()) {
        m_log << MSG::ERROR << "could not retreive the CondAttrListCollection from DB folder " << m_deadPanelFolder << endmsg;
        return sc;
    }

    else if (m_debug)
        m_log << MSG::DEBUG << " CondAttrListCollection from DB folder have been obtained with size " << atrc->size() << endmsg;

    CondAttrListCollection::const_iterator itr;
    for (itr = atrc->begin(); itr != atrc->end(); ++itr) {
        // itr=atrc->chanAttrListPair(chanNum);
        const coral::AttributeList& atr = itr->second;

        std::string panel_dead;
        std::string panel_reason_dead;

        //    if(atr.size()==1){
        if (atr.size()) {
            panel_dead = *(static_cast<const std::string*>((atr["RpcDeadROPanelIds"]).addressOfData()));
            panel_reason_dead = *(static_cast<const std::string*>((atr["RpcDeadROPanelReasons"]).addressOfData()));

            if (m_debug) m_log << MSG::DEBUG << "panel_dead " << panel_dead << endmsg;
            if (m_debug) m_log << MSG::DEBUG << "panel_reason " << panel_reason_dead << endmsg;

            char delimiter = ',';
            const auto info_panel = MuonCalib::MdtStringUtils::tokenize(panel_dead, delimiter);

            Identifier PanelId;

            for (unsigned int i = 0; i < info_panel.size(); i++) {
                const auto& ch_tmp = info_panel[i];
                if (m_debug) m_log << MSG::DEBUG << " info_panel " << ch_tmp << " " << MuonCalib::MdtStringUtils::atoi(ch_tmp) << endmsg;

                PanelId = MuonCalib::MdtStringUtils::atoi(ch_tmp);

                if (PanelId.get_compact()) {
                    if (m_debug) m_log << MSG::DEBUG << "DEADPANEL " << m_rpcIdHelper->show_to_string(PanelId) << endmsg;
                    Identifier atlasId = m_rpcIdHelper->panelID(PanelId);
                    if (atlasId != 0) m_cachedDeadPanelId.push_back(PanelId);
                    if (m_debug) m_log << MSG::DEBUG << "push-back" << endmsg;
                }
            }
        }
    }

    return StatusCode::SUCCESS;
}
