/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonTrackScoringTool.h"

#include <cassert>

#include "TrkDetElementBase/TrkDetElementBase.h"
#include "TrkEventPrimitives/FitQuality.h"
#include "TrkRIO_OnTrack/RIO_OnTrack.h"
#include "TrkTrack/Track.h"
#include "TrkTrackSummary/TrackSummary.h"

namespace Muon {

    MuonTrackScoringTool::MuonTrackScoringTool(const std::string& t, const std::string& n, const IInterface* p) :
        AthAlgTool(t, n, p), m_summaryTypeScore(Trk::numberOfTrackSummaryTypes) {
        declareInterface<Trk::ITrackScoringTool>(this);

        // set some test values
        m_summaryTypeScore[Trk::numberOfPixelHits] = 20;
        m_summaryTypeScore[Trk::numberOfPixelSharedHits] = -10;  // a shared hit is only half the weight
        m_summaryTypeScore[Trk::numberOfPixelHoles] = -10;       // a hole is bad

        m_summaryTypeScore[Trk::numberOfInnermostPixelLayerHits] = 10;        // addition for being b-layer
        m_summaryTypeScore[Trk::numberOfInnermostPixelLayerSharedHits] = -5;  // a shared hit is only half the weight

        m_summaryTypeScore[Trk::numberOfGangedPixels] = -5;  // decrease for being ganged

        m_summaryTypeScore[Trk::numberOfSCTHits] = 10;        // half of a pixel, since only 1dim
        m_summaryTypeScore[Trk::numberOfSCTSharedHits] = -5;  // a shared hit is only half the weight
        m_summaryTypeScore[Trk::numberOfSCTHoles] = -5;       // a hole is bad !

        m_summaryTypeScore[Trk::numberOfTRTHits] = 2;               // 5 straws ~ 1 SCT
        m_summaryTypeScore[Trk::numberOfTRTHighThresholdHits] = 1;  // addition for being TR
        m_summaryTypeScore[Trk::numberOfOutliersOnTrack] = -2;      // an outlier might happen

        // scoring for Muons is missing
        m_summaryTypeScore[Trk::numberOfMdtHits] = 5;
        m_summaryTypeScore[Trk::numberOfTgcPhiHits] = 5;
        m_summaryTypeScore[Trk::numberOfTgcEtaHits] = 5;
        m_summaryTypeScore[Trk::numberOfCscPhiHits] = 5;
        m_summaryTypeScore[Trk::numberOfCscEtaHits] = 5;
        m_summaryTypeScore[Trk::numberOfRpcPhiHits] = 5;
        m_summaryTypeScore[Trk::numberOfRpcEtaHits] = 5;
        // New Small Wheel
        m_summaryTypeScore[Trk::numberOfStgcPhiHits] = 5;
        m_summaryTypeScore[Trk::numberOfStgcEtaHits] = 5;
        m_summaryTypeScore[Trk::numberOfMmHits] = 5;
    }

    StatusCode MuonTrackScoringTool::initialize() {
        ATH_CHECK(m_trkSummaryTool.retrieve());
        ATH_MSG_DEBUG("Retrieved tool " << m_trkSummaryTool);

        return StatusCode::SUCCESS;
    }

    Trk::TrackScore MuonTrackScoringTool::score(const Trk::Track& track, const bool /*suppressHoleSearch*/) const {
        Trk::TrackScore score;
        const Trk::TrackSummary* summary = track.trackSummary();
        if (summary) {
            score = simpleScore(track, *summary);
        } else {
            // This is potentially slow, so might need revisiting.
            std::unique_ptr<Trk::TrackSummary> tmpSummary = m_trkSummaryTool->summaryNoHoleSearch(track);
            score = simpleScore(track, *tmpSummary);
        }
        return score;
    }

    Trk::TrackScore MuonTrackScoringTool::simpleScore(const Trk::Track& track, const Trk::TrackSummary& trackSummary) const {
        // --- reject bad tracks
        if (track.fitQuality() && track.fitQuality()->numberDoF() < 0) {
            ATH_MSG_VERBOSE("numberDoF < 0, reject it");
            return Trk::TrackScore(0);
        }

        ATH_MSG_DEBUG(m_printer->print(track));

        // --- now start scoring
        Trk::TrackScore score(200);  // score of 100 per track

        // --- prob(chi2,NDF), protect for chi2<0
        if (track.fitQuality() != nullptr && track.fitQuality()->chiSquared() > 0 && track.fitQuality()->numberDoF() > 0) {
            score += 5 * track.fitQuality()->numberDoF() - track.fitQuality()->chiSquared();
        }

        // --- summary score analysis
        for (int i = 0; i < Trk::numberOfTrackSummaryTypes; ++i) {
            int value = trackSummary.get(static_cast<Trk::SummaryType>(i));
            // value is -1 if undefined.
            if (value > 0) {
                score += m_summaryTypeScore[i] * value;
                ATH_MSG_VERBOSE("\tType [" << i << "], value \t= " << value << "], score \t=" << score);
            }
        }
        if (score == 0)
            score =
                0.000001;  // since 0 is the bad track score; you'd have to get very unlucky to get a score of exactly 0 but it can happen
        ATH_MSG_DEBUG(" Track Score " << score);

        return score;
    }
}  // namespace Muon
