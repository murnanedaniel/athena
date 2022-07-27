/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonCombinedInDetExtensionAlg.h"

#include "MuonCombinedEvent/MuonCandidateCollection.h"

MuonCombinedInDetExtensionAlg::MuonCombinedInDetExtensionAlg(const std::string& name, ISvcLocator* pSvcLocator) :
    AthReentrantAlgorithm(name, pSvcLocator) {}

StatusCode MuonCombinedInDetExtensionAlg::initialize() {
    ATH_MSG_VERBOSE(" usePRDs = " << m_usePRDs);
    ATH_CHECK(m_muonCombinedInDetExtensionTools.retrieve());
    ATH_CHECK(m_indetCandidateCollectionName.initialize());
    ATH_CHECK(m_MDT_ContainerName.initialize(m_usePRDs));
    ATH_CHECK(m_RPC_ContainerName.initialize(m_usePRDs));
    ATH_CHECK(m_TGC_ContainerName.initialize(m_usePRDs));
    ATH_CHECK(m_CSC_ContainerName.initialize(m_usePRDs && m_hasCSC));
    ATH_CHECK(m_sTGC_ContainerName.initialize(m_usePRDs && m_hasSTGC));
    ATH_CHECK(m_MM_ContainerName.initialize(m_usePRDs && m_hasMM));
    ATH_CHECK(m_tagMap.initialize());
    ATH_CHECK(m_combTracks.initialize(!m_combTracks.key().empty()));
    ATH_CHECK(m_METracks.initialize(!m_METracks.key().empty()));
    ATH_CHECK(m_segments.initialize(!m_segments.key().empty()));

    return StatusCode::SUCCESS;
}

StatusCode MuonCombinedInDetExtensionAlg::execute(const EventContext& ctx) const {
    SG::ReadHandle<InDetCandidateCollection> indetCandidateCollection(m_indetCandidateCollectionName, ctx);
    if (!indetCandidateCollection.isValid()) {
        ATH_MSG_ERROR("Could not read " << m_indetCandidateCollectionName);
        return StatusCode::FAILURE;
    }

    ATH_MSG_VERBOSE("Loaded InDetCandidateCollection " << m_indetCandidateCollectionName << " with  " << indetCandidateCollection->size()
                                                       << " elements.");
    if (msgLvl(MSG::VERBOSE)) {
        for (const MuonCombined::InDetCandidate* candidate : *indetCandidateCollection)
            msg(MSG::VERBOSE) << candidate->toString() << endmsg;
    }

    SG::WriteHandle<MuonCombined::InDetCandidateToTagMap> tagMap(m_tagMap, ctx);
    ATH_CHECK(tagMap.record(std::make_unique<MuonCombined::InDetCandidateToTagMap>()));

    TrackCollection* combTracks = nullptr;
    TrackCollection* meTracks = nullptr;
    Trk::SegmentCollection* segments = nullptr;

    if (!m_combTracks.key().empty()) {
        SG::WriteHandle<TrackCollection> wh_combTracks(m_combTracks, ctx);
        ATH_CHECK(wh_combTracks.record(std::make_unique<TrackCollection>()));
        combTracks = wh_combTracks.ptr();
    }
    if (!m_METracks.key().empty()) {
        SG::WriteHandle<TrackCollection> wh_meTracks(m_METracks, ctx);
        ATH_CHECK(wh_meTracks.record(std::make_unique<TrackCollection>()));
        meTracks = wh_meTracks.ptr();
    }

    if (!m_segments.key().empty()) {
        SG::WriteHandle<Trk::SegmentCollection> wh_segs(m_segments, ctx);
        ATH_CHECK(wh_segs.record(std::make_unique<Trk::SegmentCollection>()));
        segments = wh_segs.ptr();
    }

    if (m_usePRDs) {
        MuonCombined::IMuonCombinedInDetExtensionTool::MuonPrdData prdData;
        SG::ReadHandle<Muon::MdtPrepDataContainer> mdtPRDContainer(m_MDT_ContainerName, ctx);
        prdData.mdtPrds = mdtPRDContainer.cptr();
        if (m_hasCSC) {
            SG::ReadHandle<Muon::CscPrepDataContainer> cscPRDContainer(m_CSC_ContainerName, ctx);
            prdData.cscPrds = cscPRDContainer.cptr();
        }
        if (m_hasSTGC && m_hasMM) {
            SG::ReadHandle<Muon::sTgcPrepDataContainer> stgcPRDContainer(m_sTGC_ContainerName, ctx);
            SG::ReadHandle<Muon::MMPrepDataContainer> mmPRDContainer(m_MM_ContainerName, ctx);
            prdData.stgcPrds = stgcPRDContainer.cptr();
            prdData.mmPrds = mmPRDContainer.cptr();
        }
        SG::ReadHandle<Muon::RpcPrepDataContainer> rpcPRDContainer(m_RPC_ContainerName, ctx);
        prdData.rpcPrds = rpcPRDContainer.cptr();
        SG::ReadHandle<Muon::TgcPrepDataContainer> tgcPRDContainer(m_TGC_ContainerName, ctx);
        prdData.tgcPrds = tgcPRDContainer.cptr();
        for (const auto& tool : m_muonCombinedInDetExtensionTools) {
            tool->extendWithPRDs(*indetCandidateCollection, tagMap.ptr(), prdData, combTracks, meTracks, segments, ctx);
        }
    } else {
        for (const auto& tool : m_muonCombinedInDetExtensionTools) {
            tool->extend(*indetCandidateCollection, tagMap.ptr(), combTracks, meTracks, segments, ctx);
        }
    }
    return StatusCode::SUCCESS;
}
