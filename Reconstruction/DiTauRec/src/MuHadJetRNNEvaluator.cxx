/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "DiTauRec/MuHadJetRNNEvaluator.h"
#include "xAODTau/TauTrack.h"

#include <algorithm>

#include "PathResolver/PathResolver.h"
#include "tauRecTools/TauJetRNN.h"


MuHadJetRNNEvaluator::MuHadJetRNNEvaluator(const std::string &name)
    : TauRecToolBase(name), m_net_1p(nullptr), m_net_3p(nullptr) {
    // Network weight files for 1- and 3-prong taus
    // If the filename is an empty string a default value is decorated
    declareProperty("NetworkFile1P", m_weightfile_1p = "");
    declareProperty("NetworkFile3P", m_weightfile_3p = "");
    declareProperty("OutputVarname", m_output_varname = "RNNJetScore");
    declareProperty("MinChargedTracks", m_min_charged_tracks = 1);
    declareProperty("MaxTracks", m_max_tracks = 10);
    declareProperty("MaxClusters", m_max_clusters = 6);
    declareProperty("MaxClusterDR", m_max_cluster_dr = 1.0f);

    // Naming conventions for the network weight files:
    declareProperty("InputLayerScalar", m_input_layer_scalar = "scalar");
    declareProperty("InputLayerTracks", m_input_layer_tracks = "tracks");
    declareProperty("InputLayerClusters", m_input_layer_clusters = "clusters");
    declareProperty("OutputLayer", m_output_layer = "rnnid_output");
    declareProperty("OutputNode", m_output_node = "sig_prob");
    declareProperty("TrkClassifyDone", m_classifierDone = false ) ;

}

MuHadJetRNNEvaluator::~MuHadJetRNNEvaluator() {}

StatusCode MuHadJetRNNEvaluator::initialize() {
    ATH_MSG_INFO("Initializing MuHadJetRNNEvaluator with output name " << m_output_varname );

    // Use PathResolver to search for the weight files
    if (!m_weightfile_1p.empty()) {
        auto weightfile_1p = find_file(m_weightfile_1p);
        if (weightfile_1p.empty()) {
            ATH_MSG_ERROR("Could not find network weights: "
                          << m_weightfile_1p);
            return StatusCode::FAILURE;
        } else {
            ATH_MSG_INFO("Using network config [1-prong]: " << weightfile_1p);
        }
        m_weightfile_1p = weightfile_1p;
    }

    if (!m_weightfile_3p.empty()) {
        auto weightfile_3p = find_file(m_weightfile_3p);
        if (weightfile_3p.empty()) {
            ATH_MSG_ERROR("Could not find network weights: "
                          << m_weightfile_3p);
            return StatusCode::FAILURE;
        } else {
            ATH_MSG_INFO("Using network config [3-prong]: " << weightfile_3p);
        }
        m_weightfile_3p = weightfile_3p;
    }

    // Set the layer and node names in the weight file
    TauJetRNN::Config config;
    config.input_layer_scalar = m_input_layer_scalar;
    config.input_layer_tracks = m_input_layer_tracks;
    config.input_layer_clusters = m_input_layer_clusters;
    config.output_layer = m_output_layer;
    config.output_node = m_output_node;

    // Load the weights and create the network
    m_net_1p = std::make_unique<TauJetRNN>(m_weightfile_1p, config);
    if (!m_net_1p) {
        ATH_MSG_WARNING("No network configured for 1-prong taus. "
                        "Decorating defaults...");
    }

    m_net_3p = std::make_unique<TauJetRNN>(m_weightfile_3p, config);
    if (!m_net_3p) {
        ATH_MSG_WARNING("No network configured for multi-prong taus. "
                        "Decorating defaults...");
    }

    if ( m_classifierDone ) m_isoTrackType = xAOD::TauJetParameters::classifiedIsolation ;
    else m_isoTrackType = xAOD::TauJetParameters::modifiedIsolationTrack ;

    return StatusCode::SUCCESS;
}

StatusCode MuHadJetRNNEvaluator::execute(xAOD::TauJet &tau) {
    // Output variable accessor
    const SG::AuxElement::Accessor<float> output(m_output_varname);

    // Set default score and overwrite later
    output(tau) = -1111.0f;

    // Only apply to taus exceeding the configured minimum number of tracks
    const auto nTracksCharged = tau.nTracksCharged();

    if (nTracksCharged < m_min_charged_tracks) {
        return StatusCode::SUCCESS;
    }

    // Get input objects
    std::vector<const xAOD::TauTrack *> tracks;
    ATH_CHECK(get_tracks(tau, tracks));
    std::vector<const xAOD::CaloCluster *> clusters;
    ATH_CHECK(get_clusters(tau, clusters));

    // Evaluate networks
    if (nTracksCharged <= 1 && m_net_1p) {
        output(tau) = m_net_1p->compute(tau, tracks, clusters);
    }

    if (nTracksCharged > 1 && m_net_3p) {
        output(tau) = m_net_3p->compute(tau, tracks, clusters);
    }

    return StatusCode::SUCCESS;
}

TauJetRNN *MuHadJetRNNEvaluator::get_rnn_1p() {
    return m_net_1p.get();
}

TauJetRNN *MuHadJetRNNEvaluator::get_rnn_3p() {
    return m_net_3p.get();
}

StatusCode MuHadJetRNNEvaluator::get_tracks(
    const xAOD::TauJet &tau, std::vector<const xAOD::TauTrack *> &out) 
{
    // basing on the assumption that muon tracks have been removed in uptream algos.
    auto tracks = tau.allTracks();

    std::vector<const xAOD::TauTrack *> tauTracks_noMuon ;
    for ( auto  ttItr :  tracks )
    {

      bool coreTrack = ttItr->flagWithMask( ( 1 << xAOD::TauJetParameters::TauTrackFlag::classifiedCharged ) ) ;
      bool modifIso = ttItr->flagWithMask(  ( 1 << ( static_cast<xAOD::TauJetParameters::TauTrackFlag>( m_isoTrackType ) ) ) ) ;

      if ( ! ( coreTrack || modifIso ) ) continue ;

      tauTracks_noMuon.push_back( ttItr ) ;
    }

    // Sort by descending pt
    auto cmp_pt = [](const xAOD::TauTrack *lhs, const xAOD::TauTrack *rhs) {
        return lhs->pt() > rhs->pt();
    };

    std::sort( tauTracks_noMuon.begin(), tauTracks_noMuon.end(), cmp_pt);

    // Truncate tracks
    if ( tauTracks_noMuon.size() > m_max_tracks) { tauTracks_noMuon.resize(m_max_tracks); }
    out = std::move( tauTracks_noMuon );

    return StatusCode::SUCCESS ;
}

StatusCode MuHadJetRNNEvaluator::get_clusters(
    const xAOD::TauJet &tau, std::vector<const xAOD::CaloCluster *> &out) 
{
    std::vector<const xAOD::CaloCluster *> clusters;

    const xAOD::Jet *jet_seed = *tau.jetLink();
    if (!jet_seed) {
        ATH_MSG_ERROR("Tau jet link is invalid.");
        return StatusCode::FAILURE;
    }

    static const SG::AuxElement::ConstAccessor< std::vector< double > > accMuonCluster( "overlapMuonCluster" );  
    std::vector< double > muCluster_v4 = accMuonCluster( tau ) ;
    TLorentzVector muCluster ;
    muCluster.SetPtEtaPhiE( muCluster_v4[0], muCluster_v4[1], muCluster_v4[2], muCluster_v4[3] ) ;

    for (const auto jc : jet_seed->getConstituents()) 
    {

        auto cl = dynamic_cast<const xAOD::CaloCluster *>(jc->rawConstituent());
        if (!cl) {
            ATH_MSG_ERROR("Calorimeter cluster is invalid.");
            return StatusCode::FAILURE;
        }

        TLorentzVector clusP4 = cl->p4() ;
        // Select clusters in cone centered on the tau detector axis
        const auto lc_p4 = tau.p4(xAOD::TauJetParameters::DetectorAxis);
        if ( lc_p4.DeltaR( clusP4 ) >= m_max_cluster_dr) continue ;

        if (     muCluster.Pt() > 0.  
             &&  muCluster.DeltaR( clusP4 ) < 0.05 
             && std::abs( muCluster.Pt() - clusP4.Pt() )/clusP4.Pt() < 0.2 ) 
        {
          ATH_MSG_DEBUG(  " overlapping muon Cluster found in MuHadJetRNNEvaluator::get_clusters" ) ;
          continue ;
        }

        clusters.push_back(cl);
    }

    // Sort by descending et
    auto et_cmp = [](const xAOD::CaloCluster *lhs,
                     const xAOD::CaloCluster *rhs) {
        return lhs->et() > rhs->et();
    };
    std::sort(clusters.begin(), clusters.end(), et_cmp);

    // Truncate clusters
    if (clusters.size() > m_max_clusters) {
        clusters.resize(m_max_clusters);
    }
    out = std::move(clusters);

    return StatusCode::SUCCESS ;
}

