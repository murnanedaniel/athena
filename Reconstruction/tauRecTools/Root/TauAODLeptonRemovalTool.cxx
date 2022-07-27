/*
    Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "tauRecTools/TauAODLeptonRemovalTool.h"

TauAODLeptonRemovalTool::TauAODLeptonRemovalTool(const std::string& name):
    TauRecToolBase(name) {
}

StatusCode TauAODLeptonRemovalTool::initialize() {
    ATH_CHECK(m_elecInputContainer.initialize());
    ATH_CHECK(m_muonInputContainer.initialize());
    m_elecWpStr = m_strElecIdWpPrefix.value() + m_strMinElecIdWp.value();
    m_muonWpUi  = m_mapMuonIdWp.at(m_strMinMuonIdWp);
    return StatusCode::SUCCESS;
}

StatusCode TauAODLeptonRemovalTool::execute(xAOD::TauJet& tau) const {
    // Read in elec and muon container
    SG::ReadHandle<xAOD::ElectronContainer> elec_input_handle(m_elecInputContainer);
    SG::ReadHandle<xAOD::MuonContainer> muon_input_handle(m_muonInputContainer);
    if (bool fail_elec = !elec_input_handle.isValid(), fail_muon = !muon_input_handle.isValid(); fail_elec || fail_muon) {
        ATH_MSG_ERROR(  (fail_elec ? "Could not retrieve Electron container with key " + elec_input_handle.key() : "") +
                        (fail_muon ? "\tCould not retrieve Muon container with key " + muon_input_handle.key() : "")
    );
        return StatusCode::FAILURE;
    }
    auto elec_container = elec_input_handle.cptr();
    auto muon_container = muon_input_handle.cptr();
    //Add the Aux element as empty vector
    const SG::AuxElement::Accessor<std::vector<ElementLink<xAOD::MuonContainer>>> acc_removed_muons("removedMuons");
    const SG::AuxElement::Accessor<std::vector<ElementLink<xAOD::ElectronContainer>>> acc_removed_elecs("removedElecs");
    acc_removed_muons(tau).clear();
    acc_removed_elecs(tau).clear();
    //get the muon and electron tracks and clusters
    auto elec_and_tracks   = decltype((getElecAndTrk)(tau, *elec_container))();
    auto elec_and_clusters = decltype((getElecAndCls)(tau, *elec_container))();
    auto muon_and_tracks   = decltype((getMuonAndTrk)(tau, *muon_container))();
    auto muon_and_clusters = decltype((getMuonAndCls)(tau, *muon_container))();
    if(m_doElecTrkRm) elec_and_tracks   = getElecAndTrk(tau, *elec_container);
    if(m_doElecClsRm) elec_and_clusters = getElecAndCls(tau, *elec_container);
    if(m_doMuonTrkRm) muon_and_tracks   = getMuonAndTrk(tau, *muon_container);
    if(m_doMuonClsRm) muon_and_clusters = getMuonAndCls(tau, *muon_container);
    // if nothing found just give up here
    if(elec_and_tracks.empty() && elec_and_clusters.empty() && muon_and_tracks.empty() && muon_and_clusters.empty()) return StatusCode::SUCCESS;
    // remove the links from the tau
    auto tau_track_links = tau.allTauTrackLinksNonConst();
    auto tau_cluster_links = tau.clusterLinks();
    auto trk_removed_muons = removeTrks(tau_track_links,    muon_and_tracks);
    auto trk_removed_elecs = removeTrks(tau_track_links,    elec_and_tracks);
    auto cls_removed_muons = removeClss(tau_cluster_links,  muon_and_clusters);
    auto cls_removed_elecs = removeClss(tau_cluster_links,  elec_and_clusters);
    tau.clearTauTrackLinks();
    tau.clearClusterLinks();
    tau.setClusterLinks(tau_cluster_links);
    tau.setAllTauTrackLinks(tau_track_links);
    //Merge the resulting vector and add them to sets
    auto removed_muons = std::move(trk_removed_muons);
    auto removed_elecs = std::move(trk_removed_elecs);
    removed_muons.insert(removed_muons.end(), cls_removed_muons.begin(), cls_removed_muons.end());
    removed_elecs.insert(removed_elecs.end(), cls_removed_elecs.begin(), cls_removed_elecs.end());
    auto removed_muons_set = std::set(removed_muons.begin(), removed_muons.end());
    auto removed_elecs_set = std::set(removed_elecs.begin(), removed_elecs.end());
    //set link to the removed lepton
    for (auto muon : removed_muons_set ){
        ElementLink<xAOD::MuonContainer> link;
        link.toContainedElement(*muon_container, muon);
        acc_removed_muons(tau).push_back(link);
    }
    for (auto elec : removed_elecs_set){
        ElementLink<xAOD::ElectronContainer> link;
        link.toContainedElement(*elec_container, elec);
        acc_removed_elecs(tau).push_back(link);
    }
    //notify the runner alg that the tau was modified
    if (!acc_removed_elecs(tau).empty() || !acc_removed_muons(tau).empty())
    {
        const SG::AuxElement::Accessor<char> acc_modified("ModifiedInAOD");
        acc_modified(tau) = static_cast<char>(true);
    }
    return StatusCode::SUCCESS;
}

//helpers
std::vector<const xAOD::CaloCluster*> TauAODLeptonRemovalTool::getOrignalTopoClusters(const xAOD::CaloCluster *cluster) const {
    static const SG::AuxElement::Accessor<std::vector<ElementLink<xAOD::CaloClusterContainer>>> acc_origClusterLinks("constituentClusterLinks");
    std::vector< const xAOD::CaloCluster* > orig_cls;
    if(acc_origClusterLinks.isAvailable(*cluster)) {
        auto links = acc_origClusterLinks(*cluster);
        for (const auto &link : links) {
            if (link.dataID() != "CaloCalTopoClusters")
                ATH_MSG_WARNING("the clusters in the lepton cannot be converted to CaloCalTopoClusters, the ID is " << link.dataID());
            if (link.isValid())
                orig_cls.push_back(*link);
        }
    }
    return orig_cls;
}

const xAOD::TrackParticle* TauAODLeptonRemovalTool::getOrignalTrackParticle(const xAOD::TrackParticle* trk) const {
    static const SG::AuxElement::Accessor<ElementLink<xAOD::TrackParticleContainer>> acc_origTracks ("originalTrackParticle");
    const xAOD::TrackParticle* orig_trk = nullptr;
    if(acc_origTracks.isAvailable(*trk)) {
        if (const auto & orig_link = acc_origTracks(*trk); orig_link.isValid()) {
            if (orig_link.dataID() != "InDetTrackParticles")
                ATH_MSG_WARNING("the tracks in the lepton cannot be converted to InDetTrackParticles, the ID is " << orig_link.dataID());
            orig_trk = *orig_link;
        }
    }
    return orig_trk;
}

std::vector<std::pair<const xAOD::TrackParticle*, const xAOD::Electron*>> TauAODLeptonRemovalTool::getElecAndTrk(const xAOD::TauJet& tau, const xAOD::ElectronContainer& elec_container) const {
    std::vector<std::pair<const xAOD::TrackParticle*, const xAOD::Electron*>> ret;
    std::for_each(elec_container.cbegin(), elec_container.cend(),
        [&](auto elec) -> void {
            if(tau.p4().DeltaR(elec->p4()) < m_lepRemovalConeSize && elec->passSelection(m_elecWpStr)) {
                auto elec_ID_tracks_links = elec->trackParticleLinks();
                for (const auto &elec_ID_tracks_link : elec_ID_tracks_links) {
                    if (elec_ID_tracks_link.isValid()) {
                        if(auto orig_ele_trk = getOrignalTrackParticle(*elec_ID_tracks_link); orig_ele_trk)
                            ret.push_back(std::make_pair(orig_ele_trk, elec));
                    }
                }
            }
        }
    );
    return ret;
}

std::vector<std::pair<const xAOD::CaloCluster*, const xAOD::Electron*>> TauAODLeptonRemovalTool::getElecAndCls(const xAOD::TauJet& tau, const xAOD::ElectronContainer& elec_container) const {
    std::vector<std::pair<const xAOD::CaloCluster*, const xAOD::Electron*>> ret;
    std::for_each(elec_container.cbegin(), elec_container.cend(),
        [&](auto elec) -> void {
            if(tau.p4().DeltaR(elec->p4()) < m_lepRemovalConeSize && elec->passSelection(m_elecWpStr)) {
                auto elec_cluster_links = elec->caloClusterLinks();
                for (const auto & elec_cluster_link : elec_cluster_links) {
                    if (elec_cluster_link.isValid()) {
                        auto orig_elec_clusters = std::move(getOrignalTopoClusters(*elec_cluster_link));
                        for (auto cluster : orig_elec_clusters){
                            ret.push_back(std::make_pair(cluster, elec));
                        }
                    }
                }
            }
        }
    );
    return ret;
}

std::vector<std::pair<const xAOD::TrackParticle*, const xAOD::Muon*>> TauAODLeptonRemovalTool::getMuonAndTrk(const xAOD::TauJet& tau, const xAOD::MuonContainer& muon_container) const {
    std::vector<std::pair<const xAOD::TrackParticle*, const xAOD::Muon*>> ret;
    std::for_each(muon_container.cbegin(), muon_container.cend(),
        [&](auto muon) -> void {
            if(tau.p4().DeltaR(muon->p4()) < m_lepRemovalConeSize && muon->quality() <= m_muonWpUi) {
                if(const auto & muon_ID_tracks_link = muon->inDetTrackParticleLink();  muon_ID_tracks_link.isValid())
                    ret.push_back(std::make_pair(std::move(*muon_ID_tracks_link), muon));
            }
        }
    );
    return ret;
}

std::vector<std::pair<const xAOD::CaloCluster*, const xAOD::Muon*>> TauAODLeptonRemovalTool::getMuonAndCls(const xAOD::TauJet& tau, const xAOD::MuonContainer& muon_container) const {
    std::vector<std::pair<const xAOD::CaloCluster*, const xAOD::Muon*>> ret;
    std::for_each(muon_container.cbegin(), muon_container.cend(),
        [&](auto muon) -> void {
            if(tau.p4().DeltaR(muon->p4()) < m_lepRemovalConeSize && muon->quality() <= m_muonWpUi) {
                if(const auto & muon_cluster_link = muon->clusterLink();  muon_cluster_link.isValid()) {
                    auto muon_cluster = std::move(*muon_cluster_link);
                    auto muon_e = muon->e();
                    auto loss_e = muon->floatParameter(xAOD::Muon::ParamEnergyLoss);
                    auto cls_e = muon_cluster->e();
                    auto loss_diff = ((cls_e - loss_e) / (cls_e + loss_e));
                    if (muon_e > cls_e && loss_diff < 0.1 && loss_diff > -0.3) {
                        auto orig_muon_clusters = std::move(getOrignalTopoClusters(muon_cluster));
                        for (auto cluster : orig_muon_clusters)
                            ret.push_back(std::make_pair(cluster, muon));
                    }
                }
            }
        }
    );
    return ret;
}

template<typename Tlep, typename Tlinks> std::vector<Tlep> TauAODLeptonRemovalTool::removeTrks(Tlinks& tau_trk_links, std::vector<std::pair<const xAOD::TrackParticle*, Tlep>>& tracks_and_leps) const {
    std::vector<Tlep> ret;
    tau_trk_links.erase(
        std::remove_if(tau_trk_links.begin(), tau_trk_links.end(),
            [&](auto tau_trk_link) -> bool {
                bool match = false;
                if(tau_trk_link.isValid()) {
                    auto tau_trk = (*tau_trk_link)->track();
                    auto where = std::find_if(tracks_and_leps.cbegin(), tracks_and_leps.cend(),
                        [&](auto track_and_lep){ return tau_trk == track_and_lep.first; });
                    if(where != tracks_and_leps.cend()) {
                        ATH_MSG_DEBUG("track with pt " << tau_trk->pt()/1000 << " GeV removed");
                        ret.push_back(where->second);
                        match = true;
                    }
                }
                return match;
            }
        ),
        tau_trk_links.end()
    );
    return ret;
}

template<typename Tlep, typename Tlinks> std::vector<Tlep> TauAODLeptonRemovalTool::removeClss(Tlinks& tau_cls_links, std::vector<std::pair<const xAOD::CaloCluster*, Tlep>>& clusters_and_leps) const {
    std::vector<Tlep> ret;
    tau_cls_links.erase(
        std::remove_if(tau_cls_links.begin(), tau_cls_links.end(),
            [&](auto tau_cls_link) -> bool {
                bool match = false;
                if(tau_cls_link.isValid()) {
                    auto tau_cls = static_cast<const xAOD::CaloCluster*>(*tau_cls_link);
                    auto where = std::find_if(clusters_and_leps.cbegin(), clusters_and_leps.cend(),
                        [&](auto cluster_and_lep){ return tau_cls == cluster_and_lep.first; });
                    if(where != clusters_and_leps.cend()) {
                        ATH_MSG_DEBUG("cluster with pt " << tau_cls->pt()/1000 << " GeV removed");
                        ret.push_back(where->second);
                        match = true;
                    }
                }
                return match;
            }
        ),
        tau_cls_links.end()
    );
    return ret;
}
