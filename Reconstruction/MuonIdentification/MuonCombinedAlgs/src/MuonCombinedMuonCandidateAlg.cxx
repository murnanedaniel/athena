/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonCombinedMuonCandidateAlg.h"

MuonCombinedMuonCandidateAlg::MuonCombinedMuonCandidateAlg(const std::string& name, ISvcLocator* pSvcLocator) :
    AthReentrantAlgorithm(name, pSvcLocator) {}

StatusCode MuonCombinedMuonCandidateAlg::initialize() {
    ATH_CHECK(m_muonCandidateTool.retrieve());
    ATH_CHECK(m_muonTrackParticleLocation.initialize());
    ATH_CHECK(m_candidateCollectionName.initialize());
    ATH_CHECK(m_msOnlyTracks.initialize());
    return StatusCode::SUCCESS;
}

StatusCode MuonCombinedMuonCandidateAlg::execute(const EventContext& ctx) const {
    // retrieve MuonSpectrometer tracks
    SG::ReadHandle<xAOD::TrackParticleContainer> muonTrackParticles(m_muonTrackParticleLocation, ctx);
    if (!muonTrackParticles.isValid()) {
        ATH_MSG_ERROR("Could not read " << m_muonTrackParticleLocation);
        return StatusCode::FAILURE;
    }
    if (!muonTrackParticles.isPresent()) {
        ATH_MSG_WARNING(m_muonTrackParticleLocation << " not present");
        return StatusCode::SUCCESS;
    }
    SG::WriteHandle<MuonCandidateCollection> muonCandidates(m_candidateCollectionName, ctx);
    ATH_CHECK(muonCandidates.record(std::make_unique<MuonCandidateCollection>()));
    SG::WriteHandle<TrackCollection> msOnlyTracks(m_msOnlyTracks, ctx);
    ATH_CHECK(msOnlyTracks.record(std::make_unique<TrackCollection>()));
    m_muonCandidateTool->create(*muonTrackParticles, *muonCandidates, *msOnlyTracks, ctx);

    return StatusCode::SUCCESS;
}
