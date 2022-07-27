/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonReadoutGeometry/MuonDetectorManager.h"

#include <TString.h>  // for Form

#include <fstream>
#include <utility>

#include "AthenaKernel/getMessageSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GeoPrimitives/GeoPrimitivesHelpers.h"
#include "MuonAlignmentData/ALinePar.h"
#include "MuonAlignmentData/BLinePar.h"
#include "MuonAlignmentData/CscInternalAlignmentPar.h"
#include "MuonReadoutGeometry/CscDetectorElement.h"
#include "MuonReadoutGeometry/CscReadoutElement.h"
#include "MuonReadoutGeometry/GlobalUtilities.h"
#include "MuonReadoutGeometry/MMReadoutElement.h"
#include "MuonReadoutGeometry/MdtDetectorElement.h"
#include "MuonReadoutGeometry/MdtReadoutElement.h"
#include "MuonReadoutGeometry/MuonStation.h"
#include "MuonReadoutGeometry/RpcDetectorElement.h"
#include "MuonReadoutGeometry/RpcReadoutElement.h"
#include "MuonReadoutGeometry/TgcDetectorElement.h"
#include "MuonReadoutGeometry/TgcReadoutElement.h"
#include "MuonReadoutGeometry/sTgcReadoutElement.h"

#ifndef SIMULATIONBASE
#include "MuonCondSvc/NSWCondUtils.h"
#endif

namespace MuonGM {

    MuonDetectorManager::MuonDetectorManager() { setName("Muon"); }

    MuonDetectorManager::~MuonDetectorManager() {
        for (unsigned int p = 0; p < m_envelope.size(); ++p) { m_envelope[p]->unref(); }
    }
    template <typename read_out, size_t N> void MuonDetectorManager::clearCache(std::array<std::unique_ptr<read_out>, N>& array) {
        for (std::unique_ptr<read_out>& ele : array) {
            if (ele) ele->clearCache();
        }
    }
    template <typename read_out, size_t N> void MuonDetectorManager::fillCache(std::array<std::unique_ptr<read_out>, N>& array) {
        for (std::unique_ptr<read_out>& ele : array) {
            if (ele) ele->fillCache();
        }
    }
    template <typename read_out, size_t N> void MuonDetectorManager::refreshCache(std::array<std::unique_ptr<read_out>, N>& array) {
        for (std::unique_ptr<read_out>& ele : array) {
            if (!ele) continue;
            ele->clearCache();
            ele->fillCache();
        }
    }
    void MuonDetectorManager::clearCache() {
        clearMdtCache();
        clearRpcCache();
        clearTgcCache();
        clearCscCache();
    }

    void MuonDetectorManager::refreshCache() {
        refreshMdtCache();
        refreshRpcCache();
        refreshTgcCache();
        refreshCscCache();
    }
    void MuonDetectorManager::refreshMdtCache() {
        // NEED to fill since FillCacheInitTime = 1 is the default now.
        refreshCache(m_mdtArray);
    }
    void MuonDetectorManager::refreshRpcCache() { refreshCache(m_rpcArray); }
    void MuonDetectorManager::refreshTgcCache() { refreshCache(m_tgcArray); }
    void MuonDetectorManager::refreshCscCache() {
        if (nCscRE()) refreshCache(m_cscArray);
    }
    void MuonDetectorManager::refreshMMCache() {
        if (nMMRE()) refreshCache(m_mmcArray);
    }
    void MuonDetectorManager::refreshsTgcCache() {
        if (nsTgcRE()) refreshCache(m_stgArray);
    }

    void MuonDetectorManager::clearMdtCache() { clearCache(m_mdtArray); }
    void MuonDetectorManager::clearRpcCache() { clearCache(m_rpcArray); }
    void MuonDetectorManager::clearTgcCache() { clearCache(m_tgcArray); }
    void MuonDetectorManager::clearCscCache() {
        if (nCscRE()) clearCache(m_cscArray);
    }
    void MuonDetectorManager::clearMMCache() {
        if (nMMRE()) clearCache(m_mmcArray);
    }
    void MuonDetectorManager::clearsTgcCache() {
        if (nsTgcRE()) clearCache(m_stgArray);
    }
    void MuonDetectorManager::fillMMCache() {
        if (nMMRE()) fillCache(m_mmcArray);
    }
    void MuonDetectorManager::fillsTgcCache() {
        if (nsTgcRE()) fillCache(m_stgArray);
    }
    void MuonDetectorManager::fillCache() {
        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
        log << MSG::INFO << "Filling cache" << endmsg;
        fillMdtCache();
        fillRpcCache();
        fillTgcCache();
        fillCscCache();
    }
    void MuonDetectorManager::fillMdtCache() { fillCache(m_mdtArray); }
    void MuonDetectorManager::fillRpcCache() { fillCache(m_rpcArray); }
    void MuonDetectorManager::fillTgcCache() { fillCache(m_tgcArray); }
    void MuonDetectorManager::fillCscCache() {
        if (nCscRE()) fillCache(m_cscArray);
    }

    unsigned int MuonDetectorManager::getNumTreeTops() const { return m_envelope.size(); }

    PVConstLink MuonDetectorManager::getTreeTop(unsigned int i) const { return m_envelope[i]; }

    void MuonDetectorManager::addTreeTop(PVLink pV) {
        pV->ref();
        m_envelope.push_back(pV);
    }

    void MuonDetectorManager::addMuonStation(MuonStation* mst) {
        std::string key = muonStationKey(mst->getStationType(), mst->getEtaIndex(), mst->getPhiIndex());
        m_MuonStationMap[key] = std::unique_ptr<MuonStation>(mst);
    }

    std::string MuonDetectorManager::muonStationKey(const std::string& stName, int statEtaIndex, int statPhiIndex) const {
        std::string key;
        if (statEtaIndex < 0)
            key = stName.substr(0, 3) + "_C_zi" + MuonGM::buildString(std::abs(statEtaIndex), 2) + "fi" +
                  MuonGM::buildString(statPhiIndex, 2);
        else
            key = stName.substr(0, 3) + "_A_zi" + MuonGM::buildString(std::abs(statEtaIndex), 2) + "fi" +
                  MuonGM::buildString(statPhiIndex, 2);
        return key;
    }

    const MuonStation* MuonDetectorManager::getMuonStation(const std::string& stName, int stEtaIndex, int stPhiIndex) const {
        std::string key = muonStationKey(stName, stEtaIndex, stPhiIndex);

        std::map<std::string, std::unique_ptr<MuonStation>>::const_iterator it = m_MuonStationMap.find(key);
        if (it != m_MuonStationMap.end())
            return (*it).second.get();
        else
            return nullptr;
    }

    MuonStation* MuonDetectorManager::getMuonStation(const std::string& stName, int stEtaIndex, int stPhiIndex) {
        std::string key = muonStationKey(stName, stEtaIndex, stPhiIndex);

        std::map<std::string, std::unique_ptr<MuonStation>>::const_iterator it = m_MuonStationMap.find(key);
        if (it != m_MuonStationMap.end())
            return (*it).second.get();
        else
            return nullptr;
    }

    void MuonDetectorManager::addRpcReadoutElement(RpcReadoutElement* x, const Identifier& id) {
        // check if RE has id as identity
        if (id != x->identify()) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addRpcReadoutElement() - Trying to add RpcReadoutElement with id %s not "
                     "matching the id assigned to the RpcReadoutElement %s",
                     __FILE__, __LINE__, m_rpcIdHelper->show_to_string(id).c_str(), m_rpcIdHelper->show_to_string(x->identify()).c_str()));
        }

        // add RE to map by RE hash
        const IdentifierHash Idhash = x->detectorElementHash();
        if (Idhash >= RpcRElMaxHash) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addRpcReadoutElement() - Trying to add RpcReadoutElement with "
                     "detector-element-hash id %d outside boundaries 0-%d",
                     __FILE__, __LINE__, (unsigned int)Idhash, RpcRElMaxHash));
        } else {
            if (m_rpcArrayByHash[Idhash]) {
                throw std::runtime_error(
                    Form("File: %s, Line: %d\nMuonDetectorManager::addRpcReadoutElement() - Trying to add RpcReadoutElement with "
                         "detector-element-hash id %d id = %s at location already taken by %s",
                         __FILE__, __LINE__, (unsigned int)Idhash, m_rpcIdHelper->show_to_string(id).c_str(),
                         m_rpcIdHelper->show_to_string(m_rpcArrayByHash[Idhash]->identify()).c_str()));
            }
            m_rpcArrayByHash[Idhash] = x;
        }
        int dbz_index{-1};
        int idx = rpcIdentToArrayIdx(id, dbz_index);
        if (m_rpcArray[idx]) {
            throw std::runtime_error(
                Form("%s:%d \nMuonDetectorManager::addRpcReadoutElement() - already stored a detector element for %s is occupied by %s ",
                     __FILE__, __LINE__, m_rpcIdHelper->show_to_string(id).c_str(),
                     m_rpcIdHelper->show_to_string(m_rpcArray[idx]->identify()).c_str()));
        }
        m_rpcArray[idx] = std::unique_ptr<RpcReadoutElement>(x);
        ++m_n_rpcRE;

        // add here the RpcDetectorElement and/or add this readoutElement to the DetectorElement
        const IdentifierHash idh = x->collectionHash();
        if (idh < RpcDetElMaxHash) {
            if (!(m_rpcDEArray[idh])) {
                m_rpcDEArray[idh] = std::make_unique<RpcDetectorElement>(nullptr, this, m_rpcIdHelper->elementID(id), idh);
                ++m_n_rpcDE;
            }
            m_rpcDEArray[idh]->addRpcReadoutElement(x, dbz_index);
        } else {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addRpcReadoutElement() - Trying to add RpcDetectorElement with "
                     "data-collection-hash id %d outside boundaries 0-%d",
                     __FILE__, __LINE__, (unsigned int)idh, RpcDetElMaxHash));
        }
    }

    const RpcReadoutElement* MuonDetectorManager::getRpcReadoutElement(const Identifier& id) const {
        int idx = rpcIdentToArrayIdx(id);
        return m_rpcArray[idx].get();
    }

    const MuonClusterReadoutElement* MuonDetectorManager::getMuonClusterReadoutElement(const Identifier& id) const {
        if (m_tgcIdHelper->is_tgc(id)) return getTgcReadoutElement(id);
        if (m_rpcIdHelper->is_rpc(id)) return getRpcReadoutElement(id);
        if (m_cscIdHelper->is_csc(id)) return getCscReadoutElement(id);
        if (m_mmIdHelper->is_mm(id)) return getMMReadoutElement(id);
        if (m_stgcIdHelper->is_stgc(id)) return getsTgcReadoutElement(id);
        return nullptr;
    }
    void MuonDetectorManager::addMMReadoutElement(MMReadoutElement* x, const Identifier& id) {
        if (id != x->identify()) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMMReadoutElement() - Trying to add MMReadoutElement with id %s not "
                     "matching the id assigned to the MMReadoutElement %s",
                     __FILE__, __LINE__, m_mmIdHelper->show_to_string(id).c_str(), m_mmIdHelper->show_to_string(x->identify()).c_str()));
        }
        const int array_idx = mmIdenToArrayIdx(id);
        if (m_mmcArray[array_idx]) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMMReadoutElement_withIdFields() - this place is taken [%d]"
                     "......... this RE cannot be added",
                     __FILE__, __LINE__, array_idx));
        }
        m_mmcArray[array_idx] = std::unique_ptr<MMReadoutElement>(x);
        ++m_n_mmcRE;
    }
    void MuonDetectorManager::addsTgcReadoutElement(sTgcReadoutElement* x, const Identifier& id) {
        if (id != x->identify()) {
            throw std::runtime_error(Form(
                "File: %s, Line: %d\nMuonDetectorManager::addsTgcReadoutElement() - Trying to add sTgcReadoutElement with id %s not "
                "matching the id assigned to the sTgcReadoutElement %s",
                __FILE__, __LINE__, m_stgcIdHelper->show_to_string(id).c_str(), m_stgcIdHelper->show_to_string(x->identify()).c_str()));
        }
        const int array_idx = stgcIdentToArrayIdx(id);
        if (m_stgArray[array_idx]) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addsTgcReadoutElement() - this place is taken [%d] "
                     "......... this RE cannot be added",
                     __FILE__, __LINE__, array_idx));
        }
        m_stgArray[array_idx] = std::unique_ptr<sTgcReadoutElement>(x);
        ++m_n_stgRE;
    }

    void MuonDetectorManager::addMdtReadoutElement(MdtReadoutElement* x, const Identifier& id) {
        if (id != x->identify()) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - Trying to add MdtReadoutElement with id %s not "
                     "matching the id assigned to the MdtReadoutElement %s",
                     __FILE__, __LINE__, m_mdtIdHelper->show_to_string(id).c_str(), m_mdtIdHelper->show_to_string(x->identify()).c_str()));
        }

        // add here the MdtReadoutElement to the array by RE hash
        // use already known RE hash
        const IdentifierHash Idhash = x->detectorElementHash();
        if (Idhash >= MdtRElMaxHash) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - Trying to add MdtReadoutElement with "
                     "detector-element-hash id %d outside boundaries 0-%d",
                     __FILE__, __LINE__, (unsigned int)Idhash, MdtRElMaxHash));
        } else {
            if (m_mdtArrayByHash[Idhash]) {
                throw std::runtime_error(
                    Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - Trying to add MdtReadoutElement with "
                         "detector-element-hash id %d id = %s at location already taken by %s",
                         __FILE__, __LINE__, (unsigned int)Idhash, m_mdtIdHelper->show_to_string(id).c_str(),
                         m_mdtIdHelper->show_to_string(m_mdtArrayByHash[Idhash]->identify()).c_str()));
            }
            m_mdtArrayByHash[Idhash] = x;
        }
        // add here the MdtDetectorElement and/or add this readoutElement to the DetectorElement
        // use already known data-collection hash
        const IdentifierHash idh = x->collectionHash();
        if (idh < MdtDetElMaxHash) {
            if (!(m_mdtDEArray[idh])) {
                m_mdtDEArray[idh] = std::make_unique<MdtDetectorElement>(nullptr, this, m_mdtIdHelper->elementID(id), idh);
                m_n_mdtDE++;
            }
            m_mdtDEArray[idh]->addMdtReadoutElement(x, m_mdtIdHelper->multilayer(id));
        } else {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - Trying to add MdtDetectorElement with "
                     "data-collection-hash id %d outside boundaries 0-%d",
                     __FILE__, __LINE__, (unsigned int)idh, MdtDetElMaxHash));
        }       
        const int arrayIdx = mdtIdentToArrayIdx(id);
        if (m_mdtArray[arrayIdx]) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - this place is taken %d current id is "
                     "%s stored id %s",
                     __FILE__, __LINE__, arrayIdx, m_mdtIdHelper->show_to_string(id).c_str(),
                     m_mdtIdHelper->show_to_string(m_mdtArray[arrayIdx]->identify()).c_str()));
        }
        m_mdtArray[arrayIdx] = std::unique_ptr<MdtReadoutElement>(x);
        ++m_n_mdtRE;
    }

    const MdtReadoutElement* MuonDetectorManager::getMdtReadoutElement(const Identifier& id) const {
        const int arrayIdx = mdtIdentToArrayIdx(id);
        return m_mdtArray[arrayIdx].get();
    }
     MdtReadoutElement* MuonDetectorManager::getMdtReadoutElement(const Identifier& id) {
        const int arrayIdx = mdtIdentToArrayIdx(id);
        return m_mdtArray[arrayIdx].get();
    }
    void MuonDetectorManager::addCscReadoutElement(CscReadoutElement* x, const Identifier& id) {
        if (id != x->identify()) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addCscReadoutElement() - Trying to add CscReadoutElement with id %s not "
                     "matching the id assigned to the CscReadoutElement %s",
                     __FILE__, __LINE__, m_cscIdHelper->show_to_string(id).c_str(), m_cscIdHelper->show_to_string(x->identify()).c_str()));
        }

        // add here RE to array by hash
        const IdentifierHash Idhash = x->detectorElementHash();
        if (Idhash >= CscRElMaxHash) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addCscReadoutElement() - Trying to add CscReadoutElement with "
                     "detector-element-hash id %d outside boundaries 0-%d",
                     __FILE__, __LINE__, (unsigned int)Idhash, CscRElMaxHash));
        } else {
            if (m_cscArrayByHash[Idhash]) {
                throw std::runtime_error(
                    Form("File: %s, Line: %d\nMuonDetectorManager::addCscReadoutElement() - Trying to add CscReadoutElement with "
                         "detector-element-hash id %d id = %s at location already taken by %s",
                         __FILE__, __LINE__, (unsigned int)Idhash, m_cscIdHelper->show_to_string(id).c_str(),
                         m_cscIdHelper->show_to_string(m_cscArrayByHash[Idhash]->identify()).c_str()));
            }
            m_cscArrayByHash[Idhash] = x;
        }

        // add here the CscDetectorElement and/or add this readoutElement to the DetectorElement
        const IdentifierHash idh = x->detectorElementHash();
        if (idh < CscDetElMaxHash) {
            if (!(m_cscDEArray[idh])) {
                m_cscDEArray[idh] = std::make_unique<CscDetectorElement>(nullptr, this, m_cscIdHelper->elementID(id), idh);
                m_n_cscDE++;
            }
            m_cscDEArray[idh]->setReadoutElement(x);
        } else {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addCscReadoutElement() - Trying to add CscDetectorElement with "
                     "data-collection-hash id %d outside boundaries 0-%d",
                     __FILE__, __LINE__, (unsigned int)idh, CscDetElMaxHash));
        }

        const int array_idx = cscIdentToArrayIdx(id);
        if (m_cscArray[array_idx]) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addCscReadoutElement() - this place is taken [%d] current id is "
                     "%s stored id %s",
                     __FILE__, __LINE__, array_idx, m_cscIdHelper->show_to_string(id).c_str(),
                     m_cscIdHelper->show_to_string(m_cscArray[array_idx]->identify()).c_str()));
        }
        m_cscArray[array_idx] = std::unique_ptr<CscReadoutElement>(x);
        ++m_n_cscRE;
    }
    const CscReadoutElement* MuonDetectorManager::getCscReadoutElement(const Identifier& id) const {
        const int array_idx = cscIdentToArrayIdx(id);
        return m_cscArray[array_idx].get();
    }
     CscReadoutElement* MuonDetectorManager::getCscReadoutElement(const Identifier& id) {
        const int array_idx = cscIdentToArrayIdx(id);
        return m_cscArray[array_idx].get();
    }
    void MuonDetectorManager::addTgcReadoutElement(TgcReadoutElement* x, const Identifier& id) {
        if (id != x->identify()) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addTgcReadoutElement() - Trying to add TgcReadoutElement with id %s not "
                     "matching the id assigned to the TgcReadoutElement %s",
                     __FILE__, __LINE__, m_tgcIdHelper->show_to_string(id).c_str(), m_tgcIdHelper->show_to_string(x->identify()).c_str()));
        }

        // add RE to array by RE hash
        const IdentifierHash Idhash = x->detectorElementHash();
        if (Idhash >= TgcRElMaxHash) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addTgcReadoutElement() - Trying to add TgcReadoutElement with "
                     "detector-element-hash id %d outside boundaries 0-%d",
                     __FILE__, __LINE__, (unsigned int)Idhash, TgcRElMaxHash));
        } else {
            if (m_tgcArrayByHash[Idhash] != nullptr) {
                throw std::runtime_error(
                    Form("File: %s, Line: %d\nMuonDetectorManager::addTgcReadoutElement() - Trying to add TgcReadoutElement with "
                         "detector-element-hash id %d id = %s at location already taken by %s",
                         __FILE__, __LINE__, (unsigned int)Idhash, m_tgcIdHelper->show_to_string(id).c_str(),
                         m_tgcIdHelper->show_to_string(m_tgcArrayByHash[Idhash]->identify()).c_str()));
            }
            m_tgcArrayByHash[Idhash] = x;
        }

        // add here the TgcDetectorElement and/or add this readoutElement to the DetectorElement
        const IdentifierHash idh = x->collectionHash();
        if (idh < TgcDetElMaxHash) {
            if (!(m_tgcDEArray[idh])) {
                m_tgcDEArray[idh] = std::make_unique<TgcDetectorElement>(nullptr, this, m_tgcIdHelper->elementID(id), idh);
                m_n_tgcDE++;
            }
            m_tgcDEArray[idh]->setReadoutElement(x);
        } else {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addTgcReadoutElement() - Trying to add TgcDetectorElement with "
                     "data-collection-hash id %d outside boundaries 0-%d",
                     __FILE__, __LINE__, (unsigned int)idh, TgcDetElMaxHash));
        }

        const int array_idx = tgcIdentToArrayIdx(id);
        if (m_tgcArray[array_idx] != nullptr) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addTgcReadoutElement() - this place is taken [%d] current id is %s "
                     "stored id %s",
                     __FILE__, __LINE__, array_idx, m_tgcIdHelper->show_to_string(id).c_str(),
                     m_tgcIdHelper->show_to_string(m_tgcArray[array_idx]->identify()).c_str()));
        }

        m_tgcArray[array_idx] = std::unique_ptr<TgcReadoutElement>(x);
        ++m_n_tgcRE;
    }

    const TgcReadoutElement* MuonDetectorManager::getTgcReadoutElement(const Identifier& id) const {
        const int array_idx = tgcIdentToArrayIdx(id);
        return m_tgcArray[array_idx].get();
    }
    TgcReadoutElement* MuonDetectorManager::getTgcReadoutElement(const Identifier& id) {
        const int array_idx = tgcIdentToArrayIdx(id);
        return m_tgcArray[array_idx].get();
    }   
    const MMReadoutElement* MuonDetectorManager::getMMReadoutElement(const Identifier& id) const {
        const int array_idx = mmIdenToArrayIdx(id);
        return m_mmcArray[array_idx].get();
    }
    const sTgcReadoutElement* MuonDetectorManager::getsTgcReadoutElement(const Identifier& id) const {
        const int array_idx = stgcIdentToArrayIdx(id);
        return m_stgArray[array_idx].get();
    }   
    int MuonDetectorManager::mmIdenToArrayIdx(const Identifier& id) const {
        return mmIdenToArrayIdx(m_mmIdHelper->isSmall(id), m_mmIdHelper->stationEta(id), m_mmIdHelper->stationPhi(id),
                                m_mmIdHelper->multilayer(id));
    }
    int MuonDetectorManager::mmIdenToArrayIdx(const int isSmall, const int stEta, const int stPhi, const int ml) const {
        const int steta_index = stEta + NMMcStEtaOffset - (stEta > 0);
        const int stphi_index = 2 * (stPhi - 1) + (isSmall == 1);
        const int ml_index = ml - 1;

        if (steta_index < 0 || steta_index >= NMMcStatEta) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::mmIdenToArrayIdx() - stEtaindex out of range %d 0-%d",
                                          __FILE__, __LINE__, steta_index, NMMcStatEta - 1));
        }
        if (stphi_index < 0 || stphi_index >= NMMcStatPhi) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::mmIdenToArrayIdx() - stPhiindex out of range %d 0-%d",
                                          __FILE__, __LINE__, stphi_index, NMMcStatPhi - 1));
        }
        if (ml_index < 0 || ml_index >= NMMcChamberLayer) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::mmIdenToArrayIdx() - ml_index out of range %d 0-%d",
                                          __FILE__, __LINE__, ml_index, NMMcChamberLayer - 1));
        }
        constexpr int C = NMMcChamberLayer;
        constexpr int BxC = NMMcStatPhi * C;
        const int array_idx = steta_index * BxC + stphi_index * C + ml_index;
        return array_idx;
    }
    int MuonDetectorManager::stgcIdentToArrayIdx(const Identifier& id) const {
        return stgcIdentToArrayIdx(m_stgcIdHelper->isSmall(id), m_stgcIdHelper->stationEta(id), m_stgcIdHelper->stationPhi(id),
                                   m_stgcIdHelper->multilayer(id));
    }
    int MuonDetectorManager::stgcIdentToArrayIdx(const int isSmall, const int stEta, const int stPhi, const int ml) const {
        /// Next the array indeces
        const int steta_index = stEta + NsTgStEtaOffset - (stEta > 0);
        const int stphi_index = 2 * (stPhi - 1) + (isSmall == 1);
        const int ml_index = ml - 1;
        if (steta_index < 0 || steta_index >= NsTgStatEta) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::stgcIdentToArrayIdx() - stEtaindex out of range %d 0-%d", __FILE__, __LINE__,
                     steta_index, NsTgStatEta - 1));
        }
        if (stphi_index < 0 || stphi_index >= NsTgStatPhi) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::stgcIdentToArrayIdx() - stPhiindex out of range %d 0-%d", __FILE__, __LINE__,
                     stphi_index, NsTgStatPhi - 1));
        }
        if (ml_index < 0 || ml_index >= NsTgChamberLayer) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::stgcIdentToArrayIdx() - ml_index out of range %d 0-%d",
                                          __FILE__, __LINE__, ml_index, NsTgChamberLayer - 1));
        }
        constexpr int C = NsTgChamberLayer;
        constexpr int BxC = NsTgStatPhi * C;
        const int array_idx = steta_index * BxC + stphi_index * C + ml_index;
        return array_idx;
    }

    int MuonDetectorManager::rpcIdentToArrayIdx(const Identifier& id) const {
        int dbl_z{-1};
        return rpcIdentToArrayIdx(id, dbl_z);
    }
    int MuonDetectorManager::rpcIdentToArrayIdx(const Identifier& id, int& dbz_index) const {
        const int stationName = m_rpcIdHelper->stationName(id);
        const int stationEta = m_rpcIdHelper->stationEta(id);
        const int doubletPhi = m_rpcIdHelper->doubletPhi(id);
        const int doubletZ = m_rpcIdHelper->doubletZ(id);
        const int doubletR = m_rpcIdHelper->doubletR(id);
        const int stname_index = rpcStationTypeIdx(stationName);
        const int steta_index = stationEta + NRpcStEtaOffset;
        const int stphi_index = m_rpcIdHelper->stationPhi(id) - 1;
        const int dbr_index = doubletR - 1;
        dbz_index = doubletZ - 1;

        // BMS 5/ |stEta|= 2 / dbR = 1 and 2 / dbZ = 3
        // BMS 6/ |stEta|= 4 / dbR = 2 / dbZ = 3
        // BMS 6/ |stEta|= 4 / dbR = 1 / dbZ = 2
        // these are the special cases where we want the rpc at doubletPhi = 2
        // to be addressed with a dbz_index=dbZ+1
        if (stname_index == RpcStatType::BMS) {
            if (std::abs(stationEta) == 2 && doubletZ == 3 && doubletPhi == 2)
                ++dbz_index;
            else if (std::abs(stationEta) == 4 && doubletR == 2 && doubletZ == 3 && doubletPhi == 2)
                ++dbz_index;
            else if (std::abs(stationEta) == 4 && doubletR == 1 && doubletZ == 2 && doubletPhi == 2)
                ++dbz_index;
        }

        if (stname_index < 0 || stname_index >= NRpcStatType) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::rpcIdentToArrayIdx() - stNameindex out of range %d 0-%d", __FILE__, __LINE__,
                     stname_index, NRpcStatType - 1));
        }
        if (steta_index < 0 || steta_index >= NRpcStatEta) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::rpcIdentToArrayIdx() - stEtaindex out of range %d 0-%d",
                                          __FILE__, __LINE__, steta_index, NRpcStatEta - 1));
        }
        if (stphi_index < 0 || stphi_index >= NRpcStatPhi) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::rpcIdentToArrayIdx() - stPhiindex out of range %d 0-%d",
                                          __FILE__, __LINE__, stphi_index, NRpcStatPhi - 1));
        }
        if (dbr_index < 0 || dbr_index >= NDoubletR) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::rpcIdentToArrayIdx() - dbr_index out of range %d 0-%d",
                                          __FILE__, __LINE__, dbr_index, NDoubletR - 1));
        }
        if (dbz_index < 0 || dbz_index >= NDoubletZ) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::rpcIdentToArrayIdx() - dbz_index out of range %d 0-%d",
                                          __FILE__, __LINE__, dbz_index, NDoubletZ - 1));
        }
        /// Unfold the array by
        /// [A][B][C][D][E]
        /// a * BxCxDxE + b * CxDxE + c*DxE +d*E +e
        constexpr int E = NDoubletZ;
        constexpr int DxE = NDoubletR * E;
        constexpr int CxDxE = NRpcStatPhi * DxE;
        constexpr int BxCxDxE = NRpcStatEta * CxDxE;
        const int arrayIdx = stname_index * BxCxDxE + steta_index * CxDxE + stphi_index * DxE + dbr_index * E + dbz_index;
        return arrayIdx;
    }
    int MuonDetectorManager::mdtIdentToArrayIdx(const Identifier& id) const {
        return mdtIdentToArrayIdx(m_mdtIdHelper->stationName(id), m_mdtIdHelper->stationEta(id), m_mdtIdHelper->stationPhi(id),
                                  m_mdtIdHelper->multilayer(id));
    }
    int MuonDetectorManager::mdtIdentToArrayIdx(const int stName, const int stEta, const int stPhi, const int ml) const {
        int stname_index = stName;
        if (stName == m_mdt_EIS_stName) {
            stname_index = NMdtStatType - 4;
        } else if (stName == m_mdt_BIM_stName) {
            stname_index = NMdtStatType - 3;
        } else if (stName == m_mdt_BME_stName) {
            stname_index = NMdtStatType - 2;
        } else if (stName == m_mdt_BMG_stName) {
            stname_index = NMdtStatType - 1;
        }
        int steta_index = stEta + NMdtStEtaOffset;
        int stphi_index = stPhi - 1;
        int ml_index = ml - 1;

        if (stname_index < 0 || stname_index >= NMdtStatType) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - stNameindex out of range %d 0-%d", __FILE__,
                     __LINE__, stname_index, NMdtStatType - 1));
        }
        if (steta_index < 0 || steta_index >= NMdtStatEta) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - stEtaindex out of range %d 0-%d", __FILE__,
                     __LINE__, steta_index, NMdtStatEta - 1));
        }
        if (stphi_index < 0 || stphi_index >= NMdtStatPhi) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - stPhiindex out of range %d 0-%d", __FILE__,
                     __LINE__, stphi_index, NMdtStatPhi - 1));
        }
        if (ml_index < 0 || ml_index >= NMdtMultilayer) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::addMdtReadoutElement() - ml_index out of range %d 0-%d",
                                          __FILE__, __LINE__, ml_index, NMdtMultilayer - 1));
        }

        /// Unfold the array by
        /// [A][B][C][D]
        /// a * BxCxD + b * CxD+ c*D +d
        constexpr int D = NMdtMultilayer;
        constexpr int CxD = NMdtStatPhi * D;
        constexpr int BxCxD = NMdtStatEta * CxD;
        const int arrayIdx = stname_index * BxCxD + steta_index * CxD + stphi_index * D + ml_index;
        return arrayIdx;
    }
    int MuonDetectorManager::tgcIdentToArrayIdx(const Identifier& id) const {
        return tgcIdentToArrayIdx(m_tgcIdHelper->stationName(id), m_tgcIdHelper->stationEta(id), m_tgcIdHelper->stationPhi(id));
    }
    int MuonDetectorManager::tgcIdentToArrayIdx(const int stationName, const int stationEta, const int stationPhi) const {
        const int stname_index = stationName + NTgcStatTypeOff;
        const int zi = stationEta;
        const int steta_index = zi + NTgcStEtaOffset - (zi > 0);
        const int stphi_index = stationPhi - 1;

        if (stname_index < 0 || stname_index >= NTgcStatType) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::tgcIdentToArrayIdx() - stNameindex out of range %d 0-%d", __FILE__, __LINE__,
                     stname_index, NTgcStatType - 1));
        }
        if (steta_index < 0 || steta_index >= NTgcStatEta) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::tgcIdentToArrayIdx() - stEtaindex out of range %d 0-%d",
                                          __FILE__, __LINE__, steta_index, NTgcStatEta - 1));
        }
        if (stphi_index < 0 || stphi_index >= NTgcStatPhi) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::tgcIdentToArrayIdx() - stPhiindex out of range %d 0-%d",
                                          __FILE__, __LINE__, stphi_index, NTgcStatPhi - 1));
        }
        /// NTgcStatType * NTgcStatEta * NTgcStatPhi
        /// Unfold the array by
        /// [A][B][C]
        /// a * BxC + b * C + c
        constexpr int C = NTgcStatPhi;
        constexpr int BxC = NTgcStatEta * C;
        return stname_index * BxC + steta_index * C + stphi_index;
    }
    int MuonDetectorManager::cscIdentToArrayIdx(const Identifier& id) const {
        return cscIdentToArrayIdx(m_cscIdHelper->stationName(id), m_cscIdHelper->stationEta(id), m_cscIdHelper->stationPhi(id),
                                  m_cscIdHelper->chamberLayer(id));
    }
    int MuonDetectorManager::cscIdentToArrayIdx(const int stName, const int stEta, const int stPhi, const int ml) const {
        const int stname_index = stName + NCscStatTypeOff;
        int steta_index = stEta + NCscStEtaOffset;
        if (steta_index == 2) steta_index = 1;
        const int stphi_index = stPhi - 1;
        const int ml_index = ml - 1;

        if (stname_index < 0 || stname_index >= NCscStatType) {
            throw std::runtime_error(
                Form("File: %s, Line: %d\nMuonDetectorManager::cscIdentToArrayIdx() - stNameindex out of range %d 0-%d", __FILE__, __LINE__,
                     stname_index, NCscStatType - 1));
        }
        if (steta_index < 0 || steta_index >= NCscStatEta) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::cscIdentToArrayIdx() - stEtaindex out of range %d 0-%d",
                                          __FILE__, __LINE__, steta_index, NCscStatEta - 1));
        }
        if (stphi_index < 0 || stphi_index >= NCscStatPhi) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::cscIdentToArrayIdx() - stPhiindex out of range %d 0-%d",
                                          __FILE__, __LINE__, stphi_index, NCscStatPhi - 1));
        }
        if (ml_index < 0 || ml_index >= NCscChamberLayer) {
            throw std::runtime_error(Form("File: %s, Line: %d\nMuonDetectorManager::cscIdentToArrayIdx() - ml_index out of range %d 0-%d",
                                          __FILE__, __LINE__, ml_index, NCscChamberLayer - 1));
        }
        constexpr int D = NCscChamberLayer;
        constexpr int CxD = NCscStatPhi * D;
        constexpr int BxCxD = NCscStatEta * CxD;
        const int array_idx = stname_index * BxCxD + steta_index * CxD + stphi_index * D + ml_index;
        return array_idx;
    }

    void MuonDetectorManager::initABlineContainers() {
        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
        m_aLineContainer.clear();
        m_bLineContainer.clear();

        if (log.level() <= MSG::DEBUG)
            log << MSG::DEBUG << "Init A/B Line Containers - pointers are <" << (uintptr_t)&m_aLineContainer << "> and <"
                << (uintptr_t)&m_bLineContainer << ">" << endmsg;

        // loop over stations to fill the A-line map at start-up
        for (auto& ist : m_MuonStationMap) {
            MuonStation* ms = ist.second.get();
            int jff = ms->getPhiIndex();
            int jzz = ms->getEtaIndex();
            std::string stType = ms->getStationType();

            ALinePar newALine;
            newALine.setAmdbId(stType, jff, jzz, 0);
            if (ms->hasALines()) {
                newALine.setParameters(ms->getALine_tras(), ms->getALine_traz(), ms->getALine_trat(), ms->getALine_rots(),
                                       ms->getALine_rotz(), ms->getALine_rott());
            } else {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "No starting A-lines for Station " << stType << " Jzz/Jff " << jzz << "/" << jff << endmsg;
                newALine.setParameters(0., 0., 0., 0., 0., 0.);
            }
            newALine.isNew(true);

            Identifier id;
            //= m_mdtIdHelper->elementID(stType, jzz, jff);
            if (m_tgcIdHelper && stType.substr(0, 1) == "T") {
                // TGC case
                int stPhi = MuonGM::stationPhiTGC(stType, jff, jzz, geometryVersion());
                int stEta = 1;            // stEta for the station is stEta for the first component chamber
                if (jzz < 0) stEta = -1;  // stEta for the station is stEta for the first component chamber
                id = m_tgcIdHelper->elementID(stType, stEta, stPhi);
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Filling A-line container with entry for key = " << m_tgcIdHelper->show_to_string(id) << endmsg;
            } else if (m_cscIdHelper && stType.substr(0, 1) == "C") {
                // CSC case
                id = m_cscIdHelper->elementID(stType, jzz, jff);
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Filling A-line container with entry for key = " << m_cscIdHelper->show_to_string(id) << endmsg;
            } else if (m_rpcIdHelper && stType.substr(0, 3) == "BML" && std::abs(jzz) == 7) {
                // RPC case
                id = m_rpcIdHelper->elementID(stType, jzz, jff, 1);
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Filling A-line container with entry for key = " << m_rpcIdHelper->show_to_string(id) << endmsg;
            } else if (m_mdtIdHelper) {
                id = m_mdtIdHelper->elementID(stType, jzz, jff);
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Filling A-line container with entry for key = " << m_mdtIdHelper->show_to_string(id) << endmsg;
            }
            m_aLineContainer.emplace(id, std::move(newALine));
            if (log.level() <= MSG::DEBUG)
                log << MSG::DEBUG << "<Filling A-line container with entry for key >" << m_mdtIdHelper->show_to_string(id) << endmsg;
        }

#ifndef SIMULATIONBASE
        if (m_stgcIdHelper && m_mmIdHelper && !(m_NSWABLineAsciiPath.empty())) {
			log << MSG::DEBUG << "Using NSW AB lines from file: " << m_NSWABLineAsciiPath << endmsg;
            ALineMapContainer writeALines;
            BLineMapContainer writeBLines;
            MuonCalib::NSWCondUtils::setNSWABLinesFromAscii(m_NSWABLineAsciiPath, writeALines, writeBLines, m_stgcIdHelper, m_mmIdHelper);
            for (auto it = writeALines.cbegin(); it != writeALines.cend(); ++it) {
                Identifier id = it->first;
                ALinePar aline = it->second;
                m_aLineContainer.emplace(id, std::move(aline));
            }

            for (auto it = writeBLines.cbegin(); it != writeBLines.cend(); ++it) {
                Identifier id = it->first;
                BLinePar bline = it->second;
                m_bLineContainer.emplace(id, std::move(bline));
            }
        }
#endif

        log << MSG::INFO << "Init A/B Line Containers - done - size is respectively " << m_aLineContainer.size() << "/"
            << m_bLineContainer.size() << endmsg;
    }

    StatusCode MuonDetectorManager::updateAlignment(const ALineMapContainer& alineData, bool isData) {
#ifdef TESTBLINES
        {
            for (auto& it : m_MuonStationMap) {
                MuonStation* station = it.second.get();
                station->setDelta_fromAline(0., 0., 0., 0., 0.,
                                            0.);  // double tras, double traz, double trat, double rots, double rotz, double rott
                if (cacheFillingFlag()) {
                    station->clearCache();
                    station->fillCache();
                } else {
                    station->refreshCache();
                }
            }
        }
#endif
        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
        if (alineData.empty()) {
            if (isData) {
                log << MSG::INFO << "Empty temporary A-line container - nothing to do here" << endmsg;
            } else {
                log << MSG::DEBUG << "Got empty A-line container (expected for MC), not applying A-lines..." << endmsg;
            }
            return StatusCode::SUCCESS;
        } else
            log << MSG::INFO << "temporary A-line container with size = " << alineData.size() << endmsg;

        // loop over the container of the updates passed by the MuonAlignmentDbTool
        unsigned int nLines = 0;
        unsigned int nUpdates = 0;
        for (const auto& [ALineId, ALine] : alineData) {
            nLines++;
            std::string stType{""};
            int jff{0}, jzz{0}, job{0};
            ALine.getAmdbId(stType, jff, jzz, job);

            if (!ALine.isNew()) {
                log << MSG::WARNING << "ALinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " is not new *** skipping" << endmsg;
                continue;
            }            
            
            //if (log.level() <= MSG::DEBUG)
                log << MSG::DEBUG << "ALinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " is new. ID = " << m_mdtIdHelper->show_to_string(ALineId) << endmsg;

            //********************
            // NSW Cases
            //********************
            
            if (stType[0] == 'M' || stType[0] == 'S') {
            
                if (!nMMRE() || !nsTgcRE()) {
                    log << MSG::WARNING << "Unable to set A-line; the manager does not contain NSW readout elements" << endmsg;
                    continue;
                }
                            
                if (!m_NSWABLineAsciiPath.empty()) {
                    log << MSG::INFO << "NSW A-lines are already set via external ascii file " << m_NSWABLineAsciiPath << endmsg;
                    continue;
                }

                if (stType[0] == 'M') {
                    // Micromegas                        
                    const int array_idx  = mmIdenToArrayIdx(ALineId);
                    MMReadoutElement* RE = m_mmcArray[array_idx].get();

                    if (!RE) {
                        log << MSG::WARNING << "AlinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " *** No MM readout element found\n"
                            << "PLEASE CHECK FOR possible MISMATCHES between alignment constants from COOL and Geometry Layout in use"
                            << endmsg;
                        return StatusCode::FAILURE;
                    }
                
                    RE->setDelta(ALine);

                } else if (stType[0] == 'S') {
                    // sTGC
                    const int array_idx    = stgcIdentToArrayIdx(ALineId);
                    sTgcReadoutElement* RE = m_stgArray[array_idx].get();

                    if (!RE) {
                        log << MSG::WARNING << "AlinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " *** No sTGC readout element found\n"
                            << "PLEASE CHECK FOR possible MISMATCHES between alignment constants from COOL and Geometry Layout in use"
                            << endmsg;
                        return StatusCode::FAILURE;
                    }
                
                    RE->setDelta(ALine);
                }

                // record this A-line in the historical A-line container
                auto [it, flag] = m_aLineContainer.insert_or_assign(ALineId, ALine);
                if (log.level() <= MSG::DEBUG) {
                    if (flag)
                        log << MSG::DEBUG << "New A-line entry for Station " << stType << " at Jzz/Jff/Job " << jzz << "/" << jff << "/" << job << endmsg;
                    else 
                        log << MSG::DEBUG << "Updating existing A-line for Station " << stType << " at Jzz/Jff/Job " << jzz << "/" << jff << "/" << job << endmsg;
                }
                
                continue;
            }
             

            //********************
            // Non-NSW Cases
            //********************

            MuonStation* thisStation = getMuonStation(stType, jzz, jff);
            if (!thisStation) {
                log << MSG::WARNING << "ALinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << "*** No MuonStation found\n"
                    << "PLEASE CHECK FOR possible MISMATCHES between alignment constants from COOL and Geometry Layout in use"
                    << endmsg;
                continue;
            }

            if (job != 0) {
                // job different than 0 (standard for TGC conditions for Sept 2010 repro.)
                if (stType.substr(0, 1) == "T") {
                    if (log.level() <= MSG::DEBUG)
                        log << MSG::DEBUG << "ALinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                            << " has JOB not 0 - this is expected for TGC" << endmsg;
                } else {
                    log << MSG::WARNING << "ALinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                        << " has JOB not 0 - this is NOT EXPECTED yet for non TGC chambers - skipping this A-line" << endmsg;
                    continue;
                }
            }

            // record this A-line in the historical A-line container
            auto [it, flag] = m_aLineContainer.insert_or_assign(ALineId, ALine);
            ALinePar& newALine = it->second;
            if (flag) {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "               New entry in A-line container for Station " << stType << " at Jzz/Jff " << jzz
                        << "/" << jff << " --- in the container with key " << m_mdtIdHelper->show_to_string(ALineId) << endmsg;
            } else {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Updating extisting entry in A-line container for Station " << stType << " at Jzz/Jff " << jzz
                        << "/" << jff << endmsg;
            }

            if (job == 0) {
                float s, z, t, ths, thz, tht;
                newALine.getParameters(s, z, t, ths, thz, tht);
                if (m_controlAlines % 10 == 0) tht = 0.;
                if (int(m_controlAlines / 10) % 10 == 0) thz = 0.;
                if (int(m_controlAlines / 100) % 10 == 0) ths = 0.;
                if (int(m_controlAlines / 1000) % 10 == 0) t = 0.;
                if (int(m_controlAlines / 10000) % 10 == 0) z = 0.;
                if (int(m_controlAlines / 100000) % 10 == 0) s = 0.;
                if (m_controlAlines != 111111) newALine.setParameters(s, z, t, ths, thz, tht);
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Setting delta transform for Station " << stType << " " << jzz << " " << jff << " "
                        << " params are = " << s << " " << z << " " << t << " " << ths << " " << thz << " " << tht << endmsg;
                thisStation->setDelta_fromAline(s, z, t, ths, thz, tht);
#ifdef TESTBLINES
                newALine.setParameters(0., 0., 0., 0., 0., 0.);
                thisStation->setDelta_fromAline(0., 0., 0., 0., 0., 0.);
#endif
                if (cacheFillingFlag()) {
                    thisStation->clearCache();
                    thisStation->fillCache();
                } else {
                    thisStation->refreshCache();
                }
            } else {
                // job different than 0 (standard for TGC conditions for Sept 2010 repro.)
                float s, z, t, ths, thz, tht;
                newALine.getParameters(s, z, t, ths, thz, tht);
                if (m_controlAlines % 10 == 0) tht = 0.;
                if (int(m_controlAlines / 10) % 10 == 0) thz = 0.;
                if (int(m_controlAlines / 100) % 10 == 0) ths = 0.;
                if (int(m_controlAlines / 1000) % 10 == 0) t = 0.;
                if (int(m_controlAlines / 10000) % 10 == 0) z = 0.;
                if (int(m_controlAlines / 100000) % 10 == 0) s = 0.;
                if (m_controlAlines != 111111) newALine.setParameters(s, z, t, ths, thz, tht);
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Setting delta transform for component " << job << " of Station " << stType << " " << jzz << " "
                        << jff << " "
                        << " params are = " << s << " " << z << " " << t << " " << ths << " " << thz << " " << tht << endmsg;
                thisStation->setDelta_fromAline_forComp(job, s, z, t, ths, thz, tht);
                if (cacheFillingFlag()) {
                    thisStation->getMuonReadoutElement(job)->clearCache();
                    thisStation->getMuonReadoutElement(job)->fillCache();
                } else {
                    thisStation->getMuonReadoutElement(job)->refreshCache();
                }
            }
            nUpdates++;
        }
        log << MSG::INFO << "# of A-lines read from the ALineMapContainer in StoreGate is " << nLines << endmsg;
        log << MSG::INFO << "# of deltaTransforms updated according to A-lines         is " << nUpdates << endmsg;
        log << MSG::INFO << "# of entries in the A-lines historical container          is " << ALineContainer()->size() << endmsg;

        return StatusCode::SUCCESS;
    }

    StatusCode MuonDetectorManager::updateDeformations(const BLineMapContainer& blineData, bool isData) {
#ifdef TESTBLINES
        {
            for (auto& it : m_MuonStationMap) {
                MuonStation* station = it.second.get();
                station->clearBLineCache();
                BLinePar* BLine = new BLinePar();
                BLine->setParameters(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
                station->setBline(BLine);
                if (cacheFillingFlag()) station->fillBLineCache();
            }
        }
#endif

        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
        log << MSG::INFO << "In updateDeformations()" << endmsg;
        if (!applyMdtDeformations()) log << MSG::INFO << "Mdt deformations are disabled; will only apply NSW deformations" << endmsg;
        
        if (blineData.empty()) {
            if (isData) {
                log << MSG::INFO << "Empty temporary B-line container - nothing to do here" << endmsg;
            } else {
                log << MSG::DEBUG << "Got empty B-line container (expected for MC), not applying B-lines..." << endmsg;
            }
            return StatusCode::SUCCESS;
        } else
            log << MSG::INFO << "temporary B-line container with size = " << blineData.size() << endmsg;

        // loop over the container of the updates passed by the MuonAlignmentDbTool
        unsigned int nLines = 0;
        unsigned int nUpdates = 0;
        for (auto [BLineId, BLine] : blineData) {
            nLines++;
            std::string stType{""};
            int jff{0}, jzz{0}, job{0};
            BLine.getAmdbId(stType, jff, jzz, job);

            //********************
            // NSW Cases
            //********************
            
            if (stType[0] == 'M' || stType[0] == 'S') {
            
                if (!nMMRE() || !nsTgcRE()) {
                    log << MSG::WARNING << "Unable to set B-line; the manager does not contain NSW readout elements" << endmsg;
                    continue;
                }
                            
                if (!m_NSWABLineAsciiPath.empty()) {
                    log << MSG::INFO << "NSW B-lines are already set via external ascii file " << m_NSWABLineAsciiPath << endmsg;
                    continue;
                }

                if (stType[0] == 'M') {
                    // Micromegas                        
                    const int array_idx  = mmIdenToArrayIdx(BLineId);
                    MMReadoutElement* RE = m_mmcArray[array_idx].get();

                    if (!RE) {
                        log << MSG::WARNING << "BlinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " *** No MM readout element found\n"
                            << "PLEASE CHECK FOR possible MISMATCHES between alignment constants from COOL and Geometry Layout in use"
                            << endmsg;
                        return StatusCode::FAILURE;
                    }
                
                    RE->setBLinePar(BLine);

                } else if (stType[0] == 'S') {
                    // sTGC
                    const int array_idx    = stgcIdentToArrayIdx(BLineId);
                    sTgcReadoutElement* RE = m_stgArray[array_idx].get();

                    if (!RE) {
                        log << MSG::WARNING << "BlinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " *** No sTGC readout element found\n"
                            << "PLEASE CHECK FOR possible MISMATCHES between alignment constants from COOL and Geometry Layout in use"
                            << endmsg;
                        return StatusCode::FAILURE;
                    }
                
                    RE->setBLinePar(BLine);
                }

                // record this B-line in the historical B-line container
                auto [it, flag] = m_bLineContainer.insert_or_assign(BLineId, BLine);
                if (log.level() <= MSG::DEBUG) {
                    if (flag)
                        log << MSG::DEBUG << "New B-line entry for Station " << stType << " at Jzz/Jff/Job " << jzz << "/" << jff << "/" << job << endmsg;
                    else 
                        log << MSG::DEBUG << "Updating existing B-line for Station " << stType << " at Jzz/Jff/Job " << jzz << "/" << jff << "/" << job << endmsg;
                }
                
                continue;
            }
            
            //********************
            // MDT Cases
            //********************    

            if (!applyMdtDeformations()) continue; // nothing to more to do
        
#ifdef TESTBLINES
            BLine.setParameters(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
#endif
            if (mdtDeformationFlag() > 999999) {
                // first reset everything
                BLine.setParameters(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
                // now apply user choice
                int choice = mdtDeformationFlag();
                if (int(choice % 10) > 0)
                    BLine.setParameters(0., 0., 0., BLine.sp(), BLine.sn(), BLine.tw(), 0., 0., BLine.eg(), BLine.ep(), 100.);
                if (int(choice % 100) > 9)
                    BLine.setParameters(0., 0., 0., BLine.sp(), BLine.sn(), BLine.tw(), 0., 0., BLine.eg(), 100., BLine.en());
                if (int(choice % 1000) > 99)
                    BLine.setParameters(0., 0., 0., BLine.sp(), BLine.sn(), BLine.tw(), 0., 0., 100., BLine.ep(), BLine.en());
                if (int(choice % 10000) > 999)
                    BLine.setParameters(0., 0., 0., BLine.sp(), BLine.sn(), 100., 0., 0., BLine.eg(), BLine.ep(), BLine.en());
                if (int(choice % 100000) > 9999)
                    BLine.setParameters(0., 0., 0., BLine.sp(), 100., BLine.tw(), 0., 0., BLine.eg(), BLine.ep(), BLine.en());
                if (int(choice % 1000000) > 99999)
                    BLine.setParameters(0., 0., 0., 100., BLine.sn(), BLine.tw(), 0., 0., BLine.eg(), BLine.ep(), BLine.en());
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Testing B-lines: control flag " << choice << " hard coding Bline = ( bz=" << BLine.bz()
                        << " bp=" << BLine.bp() << " bn=" << BLine.bn() << " sp=" << BLine.sp() << " sn=" << BLine.sn()
                        << " tw=" << BLine.tw() << " pg=" << BLine.pg() << " tr=" << BLine.tr() << " eg=" << BLine.eg()
                        << " ep=" << BLine.ep() << " en=" << BLine.en() << ")" << endmsg;
            }

            if (stType.substr(0, 1) == "T" || stType.substr(0, 1) == "C" || (stType.substr(0, 3) == "BML" && std::abs(jzz) == 7)) {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "BLinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                        << " is not a MDT station - skipping" << endmsg;
                continue;
            }
            if (mdtDeformationFlag() == 2 &&
                (stType.substr(0, 3) == "BEE" || stType.substr(0, 1) == "E"))  // MDT deformations are requested for Barrel(ASAP) only !!!!
            {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << " mdtDeformationFlag()==" << mdtDeformationFlag() << " stName = " << stType.substr(0, 3)
                        << " barrel / ec initial = " << stType.substr(0, 1) << " 	 skipping this b-line" << endmsg;
                continue;  // MDT deformations are requested for Barrel(ASAP) only !!!!
            }
            if (mdtDeformationFlag() == 3 && (stType.substr(0, 3) != "BEE" && stType.substr(0, 1) == "B")) {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << " mdtDeformationFlag()==" << mdtDeformationFlag() << " stName = " << stType.substr(0, 3)
                        << " barrel / ec initial = " << stType.substr(0, 1) << " 	 skipping this b-line" << endmsg;
                continue;  // MDT deformations are requested for Endcap(ARAMYS) only !!!!
            }
            if (mdtDeformationFlag() == 0) {
                if (log.level() <= MSG::DEBUG) log << MSG::DEBUG << " mdtDeformationFlag()==0 skipping this b-line" << endmsg;
                continue;  // should never happen...
            }
            if (!BLine.isNew()) {
                log << MSG::WARNING << "BLinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                    << " is not new *** skipping" << endmsg;
                continue;
            }
            if (log.level() <= MSG::DEBUG)
                log << MSG::DEBUG << "BLinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                    << " is new ID = " << m_mdtIdHelper->show_to_string(BLineId) << endmsg;
            if (job == 0) {
                MuonStation* thisStation = getMuonStation(stType, jzz, jff);
                if (!thisStation) {
                    log << MSG::WARNING << "BLinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                        << " *** No MuonStation found \n PLEASE CHECK FOR possible MISMATCHES between alignment constants from COOL and "
                           "Geometry Layout in use"
                        << endmsg;
                    continue;
                }

                // record this B-line in the historical B-line container
                auto [it, flag] = m_bLineContainer.insert_or_assign(BLineId, BLine);
                if (flag) {
                    if (log.level() <= MSG::DEBUG)
                        log << MSG::DEBUG << "               New entry in B-line container for Station " << stType << " at Jzz/Jff " << jzz
                            << "/" << jff << " --- in the container with key " << m_mdtIdHelper->show_to_string(BLineId) << endmsg;
                } else {
                    if (log.level() <= MSG::DEBUG)
                        log << MSG::DEBUG << "Updating existing entry in B-line container for Station " << stType << " at Jzz/Jff " << jzz
                            << "/" << jff << endmsg;
                }

                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Setting deformation parameters for Station " << stType << " " << jzz << " " << jff << " "
                        << endmsg;
                thisStation->clearBLineCache();
                thisStation->setBline(&it->second);
                if (cacheFillingFlag()) thisStation->fillBLineCache();
                nUpdates++;
            } else {
                log << MSG::WARNING << "BLinePar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " has JOB not 0 "
                    << endmsg;
                return StatusCode::FAILURE;
            }
        }
        log << MSG::INFO << "# of B-lines read from the ALineMapContainer in StoreGate   is " << nLines << endmsg;
        log << MSG::INFO << "# of deform-Transforms updated according to B-lines         is " << nUpdates << endmsg;
        log << MSG::INFO << "# of entries in the B-lines historical container            is " << BLineContainer()->size() << endmsg;

        return StatusCode::SUCCESS;
    }

    void MuonDetectorManager::storeTgcReadoutParams(std::unique_ptr<const TgcReadoutParams> x) {
        m_TgcReadoutParamsVec.push_back(std::move(x));
    }

    StatusCode MuonDetectorManager::initCSCInternalAlignmentMap() {
        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");

        if (!m_useCscIlinesFromGM) {
            log << MSG::INFO << "Init of CSC I-Lines will be done via Conditions DB" << endmsg;
            m_cscALineContainer.clear();

            for (auto& ist : m_MuonStationMap) {
                MuonStation* ms = ist.second.get();
                std::string stType = ms->getStationType();
                if (stType[0] != 'C') continue;

                int jff{ms->getPhiIndex()}, jzz{ms->getEtaIndex()}, job{3};  // it's always like this for CSCs

                for (unsigned int wlay = 1; wlay < 5; ++wlay) {
                    CscInternalAlignmentPar newILine;
                    newILine.setAmdbId(stType, jff, jzz, job, wlay);
                    if (log.level() <= MSG::DEBUG)
                        log << MSG::DEBUG << "No starting I-Lines or reseting them for Station " << stType << " Jzz/Jff/Wlay " << jzz << "/"
                            << jff << "/" << wlay << endmsg;
                    // there is no way to check if the RE already has parameters set - always overwriting them.
                    newILine.setParameters(0., 0., 0., 0., 0., 0.);
                    newILine.isNew(true);
                    Identifier idp = m_cscIdHelper->parentID(ms->getMuonReadoutElement(job)->identify());
                    Identifier id = m_cscIdHelper->channelID(idp, 2, wlay, 0, 1);
                    if (log.level() <= MSG::DEBUG)
                        log << MSG::DEBUG << "<Filling I-Line container with entry for key >" << m_cscIdHelper->show_to_string(id)
                            << endmsg;
                    m_cscALineContainer.emplace(id, newILine);
                }
            }
            log << MSG::INFO << "Init I-Line Container - done - size is respectively " << m_cscALineContainer.size() << endmsg;
        }
        if (log.level() <= MSG::DEBUG)
            log << MSG::DEBUG << "Init CSC I-Line Containers - pointer is <" << (uintptr_t)&m_cscALineContainer << ">" << endmsg;

        log << MSG::INFO << "I-Line for CSC wire layers loaded (Csc Internal Alignment)" << endmsg;
        if (m_useCscIntAlign)
            log << MSG::INFO << "According to configuration they WILL be used " << endmsg;
        else
            log << MSG::INFO << "According to configuration parameters they WILL BE UPDATED FROM CONDDB " << endmsg;
        return StatusCode::SUCCESS;
    }
    StatusCode MuonDetectorManager::updateCSCInternalAlignmentMap(const CscInternalAlignmentMapContainer& ilineData) {
        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
        if (ilineData.empty()) {
            log << MSG::WARNING << "Empty temporary CSC I-line container - nothing to do here" << endmsg;
            return StatusCode::SUCCESS;
        } else
            log << MSG::INFO << "temporary CSC I-line container with size = " << ilineData.size() << endmsg;

        // loop over the container of the updates passed by the MuonAlignmentDbTool
        unsigned int nLines{0}, nUpdates{0};
        for (const auto& [ILineId, ILine] : ilineData) {
            nLines++;
            std::string stType = "";
            int jff{0}, jzz{0}, job{0}, jlay{0};
            ILine.getAmdbId(stType, jff, jzz, job, jlay);
            if (!ILine.isNew()) {
                log << MSG::WARNING << "CscInternalAlignmentPar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " "
                    << jlay << " is not new *** skipping" << endmsg;
                continue;
            }
            if (log.level() <= MSG::DEBUG)
                log << MSG::DEBUG << "CscInternalAlignmentPar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " "
                    << jlay << " is new ID = " << m_cscIdHelper->show_to_string(ILineId) << endmsg;
            if (job == 3) {
                MuonStation* thisStation = getMuonStation(stType, jzz, jff);
                if (!thisStation) {
                    log << MSG::WARNING << "CscInternalAlignmentPar with AmdbId " << stType << " " << jzz << " " << jff << " " << job << " "
                        << jlay
                        << " *** No MuonStation found \n PLEASE CHECK FOR possible MISMATCHES between alignment constants from COOL and "
                           "Geometry Layout in use"
                        << endmsg;
                    continue;
                }

                auto [it, flag] = m_cscALineContainer.insert_or_assign(ILineId, ILine);
                if (flag) {
                    if (log.level() <= MSG::DEBUG)
                        log << MSG::DEBUG << "               New entry in CSC I-line container for Station " << stType
                            << " at Jzz/Jff/Jlay " << jzz << "/" << jff << "/" << jlay << " --- in the container with key "
                            << m_cscIdHelper->show_to_string(ILineId) << endmsg;
                } else {
                    if (log.level() <= MSG::DEBUG)
                        log << MSG::DEBUG << "Updating extisting entry in CSC I-line container for Station " << stType
                            << " at Jzz/Jff/Jlay " << jzz << "/" << jff << "/" << jlay << endmsg;
                }

                CscInternalAlignmentPar& newILine = it->second;
                float tras{0.f}, traz{0.f}, trat{0.f}, rots{0.f}, rotz{0.f}, rott{0.f};
                newILine.getParameters(tras, traz, trat, rots, rotz, rott);
                int choice = CscIlinesFlag();
                if (choice % 10 == 0) tras = 0.;
                if (int(choice / 10) % 10 == 0) rotz = 0.;
                if (int(choice / 100) % 10 == 0) rots = 0.;
                if (int(choice / 1000) % 10 == 0) trat = 0.;
                if (int(choice / 10000) % 10 == 0) traz = 0.;
                if (int(choice / 100000) % 10 == 0) traz = 0.;
                if (m_controlCscIlines != 111111) newILine.setParameters(tras, traz, trat, rots, rotz, rott);
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Setting CSC I-Lines for Station " << stType << " " << jzz << " " << jff << " " << job << " "
                        << jlay << " "
                        << " params are = " << tras << " " << traz << " " << trat << " " << rots << " " << rotz << " " << rott << endmsg;
                CscReadoutElement* CscRE = dynamic_cast<CscReadoutElement*>(thisStation->getMuonReadoutElement(job));
                if (!CscRE)
                    log << MSG::ERROR << "The CSC I-lines container includes stations which are no CSCs! This is impossible." << endmsg;
                else {
                    CscRE->setCscInternalAlignmentPar(newILine);
                }
                if (cacheFillingFlag()) {
                    thisStation->clearCache();
                    thisStation->fillCache();
                } else {
                    thisStation->refreshCache();
                }
                nUpdates++;

            } else {
                log << MSG::ERROR << "job for CSC I-Lines= " << job << " is not 3 => This is not valid." << endmsg;
            }
        }
        log << MSG::INFO << "# of CSC I-lines read from the ILineMapContainer in StoreGate is " << nLines << endmsg;
        log << MSG::INFO << "# of deltaTransforms updated according to A-lines             is " << nUpdates << endmsg;
        log << MSG::INFO << "# of entries in the CSC I-lines historical container          is " << CscInternalAlignmentContainer()->size()
            << endmsg;

        return StatusCode::SUCCESS;
    }
    StatusCode MuonDetectorManager::updateMdtAsBuiltParams(const MdtAsBuiltMapContainer& asbuiltData) {
        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
        if (asbuiltData.empty()) {
            log << MSG::WARNING << "Empty temporary As-Built container - nothing to do here" << endmsg;
            return StatusCode::SUCCESS;
        } else
            log << MSG::INFO << "temporary As-Built container with size = " << asbuiltData.size() << endmsg;

        // loop over the container of the updates passed by the MuonAlignmentDbTool
        unsigned int nLines{0}, nUpdates{0};
        for (const auto& [AsBuiltId, AsBuiltPar] : asbuiltData) {
            nLines++;
            std::string stType = "";
            int jff{0}, jzz{0}, job{0};
            AsBuiltPar.getAmdbId(stType, jff, jzz, job);
            if (!AsBuiltPar.isNew()) {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "MdtAsBuiltPar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                        << " is not new *** skipping" << endmsg;
                continue;
            }

            auto [it, flag] = m_AsBuiltParamsMap.insert_or_assign(AsBuiltId, AsBuiltPar);
            if (flag) {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "New entry in AsBuilt container for Station " << stType << " at Jzz/Jff " << jzz << "/" << jff
                        << " --- in the container with key " << m_mdtIdHelper->show_to_string(AsBuiltId) << endmsg;
            } else {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Updating extisting entry in AsBuilt container for Station " << stType << " at Jzz/Jff " << jzz
                        << "/" << jff << endmsg;
            }

            if (log.level() <= MSG::DEBUG)
                log << MSG::DEBUG << "MdtAsBuiltPar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                    << " is new ID = " << m_mdtIdHelper->show_to_string(AsBuiltId) << endmsg;

            MuonStation* thisStation = getMuonStation(stType, jzz, jff);
            if (thisStation) {
                if (log.level() <= MSG::DEBUG)
                    log << MSG::DEBUG << "Setting as-built parameters for Station " << stType << " " << jzz << " " << jff << " " << endmsg;
                thisStation->clearBLineCache();
                thisStation->setMdtAsBuiltParams(&it->second);
                if (cacheFillingFlag()) thisStation->fillBLineCache();
                nUpdates++;
            } else {
                log << MSG::WARNING << "MdtAsBuiltPar with AmdbId " << stType << " " << jzz << " " << jff << " " << job
                    << " *** No MuonStation found \n PLEASE CHECK FOR possible MISMATCHES between alignment constants from COOL and "
                       "Geometry Layout in use"
                    << endmsg;
                continue;
            }
        }
        log << MSG::INFO << "# of MDT As-Built read from the MdtAsBuiltMapContainer in StoreGate is " << nLines << endmsg;
        log << MSG::INFO << "# of deltaTransforms updated according to As-Built                  is " << nUpdates << endmsg;
        log << MSG::INFO << "# of entries in the MdtAsBuilt historical container                 is " << MdtAsBuiltContainer()->size()
            << endmsg;

        return StatusCode::SUCCESS;
    }
    void MuonDetectorManager::storeCscInternalAlignmentParams(const CscInternalAlignmentPar& x) {
        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");

        std::string stName = "XXX";
        int jff{0}, jzz{0}, job{0}, wlayer{0};
        x.getAmdbId(stName, jff, jzz, job, wlayer);
        // chamberLayer is always 2 => job is always 3
        int chamberLayer = 2;
        if (job != 3)
            log << MSG::WARNING << "job = " << job << " is not 3 => chamberLayer should be 1 - not existing ! setting 2" << endmsg;
        Identifier id = m_cscIdHelper->channelID(stName, jzz, jff, chamberLayer, wlayer, 0, 1);

        m_cscALineContainer.emplace(id, x);
        if (log.level() <= MSG::DEBUG) {
            log << MSG::DEBUG << "Adding Aline for CSC wire layer: " << m_cscIdHelper->show_to_string(id) << endmsg;
            log << MSG::DEBUG << "CscInternalAlignmentMapContainer has currently size " << m_cscALineContainer.size() << endmsg;
        }
    }

    void MuonDetectorManager::storeMdtAsBuiltParams(const MdtAsBuiltPar& params) {
        MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");

        std::string stName = "XXX";
        int jff{0}, jzz{0}, job{0};
        params.getAmdbId(stName, jff, jzz, job);
        Identifier id = mdtIdHelper()->elementID(stName, jzz, jff);
        if (!id.is_valid()) {
            log << MSG::ERROR << "Invalid MDT identifiers: sta=" << stName << " eta=" << jzz << " phi=" << jff << endmsg;
            return;
        }

        if (m_AsBuiltParamsMap.insert_or_assign(id, params).second) {
            if (log.level() <= MSG::DEBUG)
                log << MSG::DEBUG << "New entry in AsBuilt container for Station " << stName << " at Jzz/Jff " << jzz << "/" << jff
                    << " --- in the container with key " << m_mdtIdHelper->show_to_string(id) << endmsg;
        } else {
            if (log.level() <= MSG::DEBUG)
                log << MSG::DEBUG << "Updating extisting entry in AsBuilt container for Station " << stName << " at Jzz/Jff " << jzz << "/"
                    << jff << endmsg;
        }

        return;
    }

    const MdtAsBuiltPar* MuonDetectorManager::getMdtAsBuiltParams(const Identifier& id) const {
        if (!MdtAsBuiltContainer()) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::DEBUG << "No Mdt AsBuilt parameter container available" << endmsg;
            return nullptr;
        }
        ciMdtAsBuiltMap iter = m_AsBuiltParamsMap.find(id);
        if (iter == m_AsBuiltParamsMap.end()) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::DEBUG << "No Mdt AsBuilt parameters for station " << id.getString()
                << " sta=" << mdtIdHelper()->stationNameString(mdtIdHelper()->stationName(id)) << " eta=" << mdtIdHelper()->stationEta(id)
                << " phi=" << mdtIdHelper()->stationPhi(id) << endmsg;
            return nullptr;
        }
        return &iter->second;
    }

    void MuonDetectorManager::setMMAsBuiltCalculator(const NswAsBuiltDbData* nswAsBuiltData) {
#ifndef SIMULATIONBASE
        m_MMAsBuiltCalculator.reset();  // unset any previous instance
        m_MMAsBuiltCalculator = std::make_unique<NswAsBuilt::StripCalculator>();
        std::string mmJson="";
        if(!nswAsBuiltData->getMmData(mmJson)){
           MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
           log << MSG::WARNING << " Cannot retrieve MM as-built conditions data from detector store!" << endmsg;
        }
        m_MMAsBuiltCalculator->parseJSON(mmJson);
#else
        // just to silence the warning about an unused parameter
        (void)nswAsBuiltData;
#endif
    }

    void MuonDetectorManager::setStgcAsBuiltCalculator(const NswAsBuiltDbData* nswAsBuiltData) {
#ifndef SIMULATIONBASE
        m_StgcAsBuiltCalculator.reset();  // unset any previous instance
        m_StgcAsBuiltCalculator = std::make_unique<NswAsBuilt::StgcStripCalculator>();
        std::string stgcJson="";
        if(!nswAsBuiltData->getSTgcData(stgcJson)){
           MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
           log << MSG::WARNING << " Cannot retrieve sTGC as-built conditions data from detector store!" << endmsg;
        }
        m_StgcAsBuiltCalculator->parseJSON(stgcJson);
#else
        // just to silence the warning about an unused parameter
        (void)nswAsBuiltData;
#endif
    }

    const MdtReadoutElement* MuonDetectorManager::getMdtReadoutElement(const IdentifierHash& id) const {
#ifndef NDEBUG
        if (id >= MdtRElMaxHash) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::WARNING << " try to getMdtReadoutElement with hashId " << (unsigned int)id << " outside range 0-"
                << MdtRElMaxHash - 1 << endmsg;
            return nullptr;
        }
#endif
        return m_mdtArrayByHash[id];
    }

    const RpcReadoutElement* MuonDetectorManager::getRpcReadoutElement(const IdentifierHash& id) const {
#ifndef NDEBUG
        if (id >= RpcRElMaxHash) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::WARNING << " try to getRpcReadoutElement with hashId " << (unsigned int)id << " outside range 0-"
                << RpcRElMaxHash - 1 << endmsg;
            return nullptr;
        }
#endif
        return m_rpcArrayByHash[id];
    }

    const TgcReadoutElement* MuonDetectorManager::getTgcReadoutElement(const IdentifierHash& id) const {
#ifndef NDEBUG
        if (id >= TgcRElMaxHash) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::WARNING << " try to getTgcReadoutElement with hashId " << (unsigned int)id << " outside range 0-"
                << TgcRElMaxHash - 1 << endmsg;
            return nullptr;
        }
#endif
        return m_tgcArrayByHash[id];
    }

    const CscReadoutElement* MuonDetectorManager::getCscReadoutElement(const IdentifierHash& id) const {
#ifndef NDEBUG
        if (id >= CscRElMaxHash) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::WARNING << " try to getCscReadoutElement with hashId " << (unsigned int)id << " outside range 0-"
                << CscRElMaxHash - 1 << endmsg;
            return nullptr;
        }
#endif
        return m_cscArrayByHash[id];
    }

    const MdtDetectorElement* MuonDetectorManager::getMdtDetectorElement(const IdentifierHash& id) const {
#ifndef NDEBUG
        if (id >= MdtDetElMaxHash) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::WARNING << " try to getMdtDetectorElement with hashId " << (unsigned int)id << " outside range 0-"
                << MdtDetElMaxHash - 1 << endmsg;
            return nullptr;
        }
#endif
        return m_mdtDEArray[id].get();
    }

    const TgcDetectorElement* MuonDetectorManager::getTgcDetectorElement(const IdentifierHash& id) const {
#ifndef NDEBUG
        if (id >= TgcDetElMaxHash) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::WARNING << " try to getTgcDetectorElement with hashId " << (unsigned int)id << " outside range 0-"
                << TgcDetElMaxHash - 1 << endmsg;
            return nullptr;
        }
#endif
        return m_tgcDEArray[id].get();
    }

    const CscDetectorElement* MuonDetectorManager::getCscDetectorElement(const IdentifierHash& id) const {
#ifndef NDEBUG
        if (id >= CscDetElMaxHash) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::WARNING << " try to getCscDetectorElement with hashId " << (unsigned int)id << " outside range 0-"
                << CscDetElMaxHash - 1 << endmsg;
            return nullptr;
        }
#endif
        return m_cscDEArray[id].get();
    }

    const RpcDetectorElement* MuonDetectorManager::getRpcDetectorElement(const IdentifierHash& id) const {
#ifndef NDEBUG
        if (id >= RpcDetElMaxHash) {
            MsgStream log(Athena::getMessageSvc(), "MGM::MuonDetectorManager");
            log << MSG::WARNING << " try to getRpcDetectorElement with hashId " << (unsigned int)id << " outside range 0-"
                << RpcDetElMaxHash - 1 << endmsg;
            return nullptr;
        }
#endif
        return m_rpcDEArray[id].get();
    }

    unsigned int MuonDetectorManager::rpcStationTypeIdx(const int stationName) const {
        std::map<int, int>::const_iterator itr = m_rpcStatToIdx.find(stationName);
        if (itr != m_rpcStatToIdx.end()) return itr->second;
        return RpcStatType::UNKNOWN;
    }

    int MuonDetectorManager::rpcStationName(const int stationIndex) const {
        std::map<int, int>::const_iterator itr = m_rpcIdxToStat.find(stationIndex);
        if (itr != m_rpcIdxToStat.end()) return itr->second;
        return -1;
    }
    void MuonDetectorManager::set_rpcIdHelper(const RpcIdHelper* idh) {
        m_rpcIdHelper = idh;
        m_rpcStatToIdx.clear();
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BML"), RpcStatType::BML));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BMS"), RpcStatType::BMS));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BOL"), RpcStatType::BOL));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BOS"), RpcStatType::BOS));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BMF"), RpcStatType::BMF));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BOF"), RpcStatType::BOF));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BOG"), RpcStatType::BOG));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BME"), RpcStatType::BME));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BIR"), RpcStatType::BIR));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BIM"), RpcStatType::BIM));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BIL"), RpcStatType::BIL));
        m_rpcStatToIdx.insert(std::pair<int, int>(m_rpcIdHelper->stationNameIndex("BIS"), RpcStatType::BIS));

        m_rpcIdxToStat.clear();
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BML, m_rpcIdHelper->stationNameIndex("BML")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BMS, m_rpcIdHelper->stationNameIndex("BMS")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BOL, m_rpcIdHelper->stationNameIndex("BOL")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BOS, m_rpcIdHelper->stationNameIndex("BOS")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BMF, m_rpcIdHelper->stationNameIndex("BMF")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BOF, m_rpcIdHelper->stationNameIndex("BOF")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BOG, m_rpcIdHelper->stationNameIndex("BOG")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BME, m_rpcIdHelper->stationNameIndex("BME")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BIR, m_rpcIdHelper->stationNameIndex("BIR")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BIM, m_rpcIdHelper->stationNameIndex("BIM")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BIL, m_rpcIdHelper->stationNameIndex("BIL")));
        m_rpcIdxToStat.insert(std::pair<int, int>(RpcStatType::BIS, m_rpcIdHelper->stationNameIndex("BIS")));
    }
    int MuonDetectorManager::mdtStationName(const int stationIndex) const {
        if (stationIndex == NMdtStatType - 4)
            return m_mdt_EIS_stName;
        else if (stationIndex == NMdtStatType - 3)
            return  m_mdt_BIM_stName;
        else if (stationIndex == NMdtStatType - 2)
            return m_mdt_BME_stName ;
        else if (stationIndex == NMdtStatType - 1)
            return m_mdt_BMG_stName;
        return stationIndex;
    }
    
    void MuonDetectorManager::setNSWABLineAsciiPath(const std::string& str) { m_NSWABLineAsciiPath = str; }
    void MuonDetectorManager::setCacheFillingFlag(int value) { m_cacheFillingFlag = value; }
    void MuonDetectorManager::setCachingFlag(int value) { m_cachingFlag = value; }
    void MuonDetectorManager::set_DBMuonVersion(const std::string& version) { m_DBMuonVersion = version; }
    void MuonDetectorManager::setGeometryVersion(const std::string& version) { m_geometryVersion = std::move(version); }
    void MuonDetectorManager::setMinimalGeoFlag(int flag) { m_minimalgeo = flag; }
    void MuonDetectorManager::setCutoutsFlag(int flag) { m_includeCutouts = flag; }
    void MuonDetectorManager::setCutoutsBogFlag(int flag) { m_includeCutoutsBog = flag; }
    void MuonDetectorManager::setGenericTgcDescriptor(const GenericTGCCache& tc) {
        m_genericTGC.frame_h = tc.frame_h;
        m_genericTGC.frame_ab = tc.frame_ab;
        m_genericTGC.nlayers = tc.nlayers;
        for (unsigned int i = 0; i < (tc.materials).size(); i++) {
            m_genericTGC.materials[i] = tc.materials[i];
            m_genericTGC.positions[i] = tc.positions[i];
            m_genericTGC.tck[i] = tc.tck[i];
        }
    }
    void MuonDetectorManager::setGenericCscDescriptor(const GenericCSCCache& cc) {
        m_genericCSC.dummy1 = cc.dummy1;
        m_genericCSC.dummy2 = cc.dummy2;
    }
    void MuonDetectorManager::setGenericMdtDescriptor(const GenericMDTCache& mc) {
        m_genericMDT.innerRadius = mc.innerRadius;
        m_genericMDT.outerRadius = mc.outerRadius;
    }
    void MuonDetectorManager::setGenericRpcDescriptor(const GenericRPCCache& rc) {
        m_genericRPC.stripSeparation = rc.stripSeparation;
        m_genericRPC.stripPanelThickness = rc.stripPanelThickness;
        m_genericRPC.rpcLayerThickness = rc.rpcLayerThickness;
        m_genericRPC.centralSupPanelThickness = rc.centralSupPanelThickness;
        m_genericRPC.GasGapThickness = rc.GasGapThickness;
        m_genericRPC.frontendBoardWidth = rc.frontendBoardWidth;
    }
    void MuonDetectorManager::set_mdtIdHelper(const MdtIdHelper* idh) { 
        m_mdtIdHelper = idh;
        m_mdt_EIS_stName = m_mdtIdHelper->stationNameIndex("EIS");
        m_mdt_BIM_stName = m_mdtIdHelper->stationNameIndex("BIM");
        m_mdt_BME_stName = m_mdtIdHelper->stationNameIndex("BME");
        m_mdt_BMG_stName = m_mdtIdHelper->stationNameIndex("BMG");
    }
    void MuonDetectorManager::set_cscIdHelper(const CscIdHelper* idh) { m_cscIdHelper = idh; }
    void MuonDetectorManager::set_tgcIdHelper(const TgcIdHelper* idh) { m_tgcIdHelper = idh; }
    void MuonDetectorManager::set_stgcIdHelper(const sTgcIdHelper* idh) { m_stgcIdHelper = idh; }
    void MuonDetectorManager::set_mmIdHelper(const MmIdHelper* idh) { m_mmIdHelper = idh; }





}  // namespace MuonGM
