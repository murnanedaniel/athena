/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetTrackSummaryHelperTool/InDetTrackSummaryHelperTool.h"

// forward declares
#include "InDetIdentifier/PixelID.h"
#include "InDetIdentifier/SCT_ID.h"
#include "InDetIdentifier/TRT_ID.h"
#include "TrkEventUtils/PRDtoTrackMap.h"
#include "TrkRIO_OnTrack/RIO_OnTrack.h"
#include "TrkTrack/Track.h"
#include "TrkTrack/TrackStateOnSurface.h"
// normal includes
#include "Identifier/Identifier.h"
#include "InDetRIO_OnTrack/PixelClusterOnTrack.h"
#include "InDetRIO_OnTrack/SCT_ClusterOnTrack.h"
#include "InDetRIO_OnTrack/TRT_DriftCircleOnTrack.h"
#include "TrkCompetingRIOsOnTrack/CompetingRIOsOnTrack.h"
#include "TrkParameters/TrackParameters.h"

#include <cassert>

//==========================================================================
InDet::InDetTrackSummaryHelperTool::InDetTrackSummaryHelperTool(
  const std::string& t,
  const std::string& n,
  const IInterface* p)
  : base_class(t, n, p)
{}

//==========================================================================

StatusCode
InDet::InDetTrackSummaryHelperTool::initialize()
{
  if (m_usePixel) {
    if (detStore()->retrieve(m_pixelId, "PixelID").isFailure()) {
      ATH_MSG_ERROR("Could not get PixelID helper !");
      return StatusCode::FAILURE;
    }
  }

  if (m_useSCT) {
    if (detStore()->retrieve(m_sctId, "SCT_ID").isFailure()) {
      ATH_MSG_ERROR("Could not get SCT_ID helper !");
      return StatusCode::FAILURE;
    }
  }

  if (m_useTRT) {
    if (detStore()->retrieve(m_trtId, "TRT_ID").isFailure()) {
      ATH_MSG_ERROR("Could not get TRT_ID helper !");
      return StatusCode::FAILURE;
    }
  }

  ATH_CHECK(
    m_assoTool.retrieve(DisableTool{ !m_doSharedHits || m_assoTool.empty() }));
  ATH_CHECK(m_holeSearchTool.retrieve(DisableTool{ m_holeSearchTool.empty() }));
  ATH_CHECK(m_TRTStrawSummaryTool.retrieve(
    DisableTool{ not m_useTRT or m_TRTStrawSummaryTool.empty() }));

  ATH_CHECK(m_clusterSplitProbContainer.initialize(
    !m_clusterSplitProbContainer.key().empty()));

  ATH_MSG_INFO("initialize() successful in " << name());

  return StatusCode::SUCCESS;
}

namespace {
bool
isShared(const Trk::PRDtoTrackMap* prd_to_track_map,
         const PublicToolHandle<Trk::IPRD_AssociationTool>& asso_tool,
         const Trk::PrepRawData& prd)
{
  if (prd_to_track_map) {
    return prd_to_track_map->isShared(prd);
  } else {
    if (!asso_tool.isEnabled()) {
      throw std::logic_error(
        "Shared hits to be computed but no PRDtoTrack provided "
        " (and no association tool configured (deprecated))");
    }
    return asso_tool->isShared(prd);
  }
}
}

//==========================================================================
void
InDet::InDetTrackSummaryHelperTool::analyse(
  const EventContext& ctx,
  const Trk::Track& track,
  const Trk::PRDtoTrackMap* prd_to_track_map,
  const Trk::RIO_OnTrack* rot,
  const Trk::TrackStateOnSurface* tsos,
  std::vector<int>& information,
  std::bitset<Trk::numberOfDetectorTypes>& hitPattern) const
{
  const Identifier& id = rot->identify();
  bool isOutlier = tsos->type(Trk::TrackStateOnSurface::Outlier);
  bool ispatterntrack = (track.info().trackFitter() == Trk::TrackInfo::Unknown);

  if (m_usePixel and m_pixelId->is_pixel(id)) {

    if (isOutlier and
        not ispatterntrack) { // ME: outliers on pattern tracks may be
                              // reintegrated by fitter, so count them as hits
      information[Trk::numberOfPixelOutliers]++;
      if (m_pixelId->layer_disk(id) == 0 and m_pixelId->is_barrel(id)) {
        information[Trk::numberOfInnermostPixelLayerOutliers]++;
      }
      if (m_pixelId->layer_disk(id) == 1 and m_pixelId->is_barrel(id)) {
        information[Trk::numberOfNextToInnermostPixelLayerOutliers]++;
      }
    } else {
      bool hitIsSplit(false);
      if (m_pixelId->is_dbm(id)) {
        int offset =
          static_cast<int>(Trk::DBM0); // get int value of first DBM layer
        offset += m_pixelId->layer_disk(id);
        hitPattern.set(offset);
        information[Trk::numberOfDBMHits]++;
      } else {
        information[Trk::numberOfPixelHits]++;
        if (m_pixelId->layer_disk(id) == 0 and m_pixelId->is_barrel(id))
          information[Trk::numberOfInnermostPixelLayerHits]++;
        if (m_pixelId->layer_disk(id) == 1 and m_pixelId->is_barrel(id))
          information[Trk::numberOfNextToInnermostPixelLayerHits]++;
        // check to see if there's an ambiguity with the ganged cluster.
        const PixelClusterOnTrack* pix = nullptr;
        if (rot->rioType(Trk::RIO_OnTrackType::PixelCluster)) {
          pix = static_cast<const PixelClusterOnTrack*>(rot);
        }
        if (not pix) {
          ATH_MSG_ERROR("Could not cast pixel RoT to PixelClusterOnTrack!");
        } else {
          const InDet::PixelCluster* pixPrd = pix->prepRawData();
          const Trk::ClusterSplitProbabilityContainer::ProbabilityInfo&
            splitProb = getClusterSplittingProbability(ctx, pixPrd);
          if (pixPrd and splitProb.isSplit()) {
            information[Trk::numberOfPixelSplitHits]++;
            hitIsSplit = true;
          }
          if (pixPrd and m_pixelId->is_barrel(id) and
              m_pixelId->layer_disk(id) == 0 and splitProb.isSplit())
            information[Trk::numberOfInnermostLayerSplitHits]++;
          if (pixPrd and m_pixelId->is_barrel(id) and
              m_pixelId->layer_disk(id) == 1 and splitProb.isSplit())
            information[Trk::numberOfNextToInnermostLayerSplitHits]++;
          if (pix->isBroadCluster())
            information[Trk::numberOfPixelSpoiltHits]++;
          if (pix->hasClusterAmbiguity()) {
            information[Trk::numberOfGangedPixels]++;
            if (pix->isFake())
              information[Trk::numberOfGangedFlaggedFakes]++;
          }
        }

        if ((m_pixelId->is_barrel(id))) {
          int offset = m_pixelId->layer_disk(id);
          if (not hitPattern.test(offset))
            information[Trk::numberOfContribPixelLayers]++;
          hitPattern.set(offset); // assumes numbered consecutively
        } else {
          int offset = static_cast<int>(
            Trk::pixelEndCap0); // get int value of first pixel endcap disc
          offset += m_pixelId->layer_disk(id);
          if (not hitPattern.test(offset))
            information[Trk::numberOfContribPixelLayers]++;
          hitPattern.set(offset); // assumes numbered consecutively
        }
      }

      if (m_doSharedHits && !isOutlier) {
        // If we are running the TIDE ambi don't count split hits as shared
        if (not(m_runningTIDE_Ambi and hitIsSplit)) {
          // used in more than one track ?
          if (isShared(prd_to_track_map, m_assoTool, *(rot->prepRawData()))) {
            ATH_MSG_DEBUG("shared Pixel hit found");
            information[Trk::numberOfPixelSharedHits]++;
            if ((m_pixelId->is_blayer(id))) {
              ATH_MSG_DEBUG("--> shared Pixel hit is in b-layer");
              information[Trk::numberOfBLayerSharedHits]++;
            }
            if ((m_pixelId->is_barrel(id) and m_pixelId->layer_disk(id) == 0)) {
              ATH_MSG_DEBUG("--> shared Pixel hit is in innermost layer");
              information[Trk::numberOfInnermostPixelLayerSharedHits]++;
            }
            if ((m_pixelId->is_barrel(id) and m_pixelId->layer_disk(id) == 1)) {
              ATH_MSG_DEBUG(
                "--> shared Pixel hit is in next to innermost layer");
              information[Trk::numberOfNextToInnermostPixelLayerSharedHits]++;
            }
          }
        }
      }
    }

  } else if (m_useSCT and m_sctId->is_sct(id)) {
    if (isOutlier and
        not ispatterntrack) { // ME: outliers on pattern tracks may be
                              // reintegrated by fitter, so count them as hits
      information[Trk::numberOfSCTOutliers]++;

    } else {
      information[Trk::numberOfSCTHits]++;

      const InDet::SCT_ClusterOnTrack* sctclus = nullptr;
      if (rot->rioType(Trk::RIO_OnTrackType::SCTCluster)) {
        sctclus = static_cast<const InDet::SCT_ClusterOnTrack*>(rot);
      }
      if (not sctclus) {
        ATH_MSG_ERROR("Could not cast SCT RoT to SCT_ClusterOnTrack!");
      } else {
        if (sctclus->isBroadCluster())
          information[Trk::numberOfSCTSpoiltHits]++;
      }

      if ((m_sctId->is_barrel(id))) {
        int offset = static_cast<int>(Trk::sctBarrel0);
        hitPattern.set(
          offset + m_sctId->layer_disk(id)); // assumes numbered consecutively
      } else {
        int offset = static_cast<int>(
          Trk::sctEndCap0); // get int value of first sct endcap disc
        hitPattern.set(
          offset + m_sctId->layer_disk(id)); // assumes numbered consecutively
      }

      if (m_doSharedHits && !isOutlier) {
        if (isShared(prd_to_track_map, m_assoTool, *(rot->prepRawData()))) {
          ATH_MSG_DEBUG("shared SCT hit found");
          information[Trk::numberOfSCTSharedHits]++;
        }
      }
    }
  } else if (m_useTRT and m_trtId->is_trt(id)) {
    bool isArgonStraw = false;
    bool isKryptonStraw = false;
    if (not m_TRTStrawSummaryTool.empty()) {
      int statusHT = m_TRTStrawSummaryTool->getStatusHT(id, ctx);
      if (statusHT == TRTCond::StrawStatus::Argon or
          statusHT == TRTCond::StrawStatus::Dead or
          statusHT == TRTCond::StrawStatus::EmulateArgon) {
        isArgonStraw = true;
      }
      if (statusHT == TRTCond::StrawStatus::Krypton or
          statusHT == TRTCond::StrawStatus::EmulateKrypton) {
        isKryptonStraw = true;
      }
    }
    if (not isArgonStraw and not isKryptonStraw) {
      information[Trk::numberOfTRTXenonHits]++;
    }

    if (isOutlier and not ispatterntrack) {
      // ME: outliers on pattern tracks may be
      // reintegrated by fitter, so count them as hits
      information[Trk::numberOfTRTOutliers]++;

      const InDet::TRT_DriftCircleOnTrack* trtDriftCircle = nullptr;
      if (rot->rioType(Trk::RIO_OnTrackType::TRT_DriftCircle)) {
        trtDriftCircle = static_cast<const InDet::TRT_DriftCircleOnTrack*>(rot);
      }
      if (not trtDriftCircle) {
        ATH_MSG_ERROR("Could not cast TRT RoT to TRT_DriftCircleOnTracknot ");
      } else {
        if (trtDriftCircle->highLevel() and not isArgonStraw and
            not isKryptonStraw)
          information[Trk::numberOfTRTHighThresholdOutliers]++;
      }
    } else {
      information[Trk::numberOfTRTHits]++;
      double error2 = rot->localCovariance()(0, 0);
      if (error2 > 1)
        information[Trk::numberOfTRTTubeHits]++;

      const InDet::TRT_DriftCircleOnTrack* trtDriftCircle = nullptr;
      if (rot->rioType(Trk::RIO_OnTrackType::TRT_DriftCircle)) {
        trtDriftCircle = static_cast<const InDet::TRT_DriftCircleOnTrack*>(rot);
      }
      if (not trtDriftCircle) {
        ATH_MSG_ERROR("Could not cast TRT RoT to TRT_DriftCircleOnTracknot ");
      } else {
        if (trtDriftCircle->highLevel()) {
          if (not isArgonStraw and not isKryptonStraw)
            information[Trk::numberOfTRTHighThresholdHits]++;
          assert(Trk::numberOfTRTHighThresholdHitsTotal < information.size());
          information[Trk::numberOfTRTHighThresholdHitsTotal]++;
        }
      }
    }

    if (m_doSharedHitsTRT && !isOutlier) {
      // used in more than one track ?
      assert(information[Trk::numberOfTRTSharedHits] >= 0);
      if (isShared(prd_to_track_map, m_assoTool, *(rot->prepRawData()))) {
        ATH_MSG_DEBUG("shared TRT hit found");
        information[Trk::numberOfTRTSharedHits]++;
      }
    }
  }
  }

void
InDet::InDetTrackSummaryHelperTool::analyse(
  const EventContext& ctx,
  const Trk::Track& track,
  const Trk::PRDtoTrackMap* prd_to_track_map,
  const Trk::CompetingRIOsOnTrack* crot,
  const Trk::TrackStateOnSurface* tsos,
  std::vector<int>& information,
  std::bitset<Trk::numberOfDetectorTypes>& hitPattern) const
{
  // re-produce prior behaviour (i.e. just take most probable ROT)
  analyse(ctx,
          track,
          prd_to_track_map,
          &crot->rioOnTrack(crot->indexOfMaxAssignProb()),
          tsos,
          information,
          hitPattern);
}

void
InDet::InDetTrackSummaryHelperTool::searchForHoles(
  const Trk::Track& track,
  std::vector<int>& information,
  const Trk::ParticleHypothesis partHyp) const
{
  ATH_MSG_DEBUG("Do hole search within HELPER, PLEASE FIX THIS AFTER 16.0.X");
  m_holeSearchTool->countHoles(track, information, partHyp);
}

void
InDet::InDetTrackSummaryHelperTool::updateSharedHitCount(
  const Trk::Track& track,
  const Trk::PRDtoTrackMap* prd_to_track_map,
  Trk::TrackSummary& summary) const
{
  // loop over track states on surface and take pixel / sct to update the shared
  // hit count
  summary.m_information[Trk::numberOfPixelSharedHits] = 0;
  summary.m_information[Trk::numberOfInnermostPixelLayerSharedHits] = 0;
  summary.m_information[Trk::numberOfNextToInnermostPixelLayerSharedHits] = 0;
  summary.m_information[Trk::numberOfSCTSharedHits] = 0;
  summary.m_information[Trk::numberOfTRTSharedHits] = 0;
  if (m_runningTIDE_Ambi) {
    summary.m_information[Trk::numberOfPixelSplitHits] = 0;
    summary.m_information[Trk::numberOfInnermostLayerSplitHits] = 0;
    summary.m_information[Trk::numberOfNextToInnermostLayerSplitHits] = 0;
  }

  const EventContext& ctx = Gaudi::Hive::currentContext();
  const DataVector<const Trk::MeasurementBase>* measurements =
    track.measurementsOnTrack();
  if (measurements) {
    for (const auto* const ms : *measurements) {
      // check if it's a rot
      const Trk::RIO_OnTrack* rot = nullptr;
      if (ms->type(Trk::MeasurementBaseType::RIO_OnTrack)) {
        rot = static_cast<const Trk::RIO_OnTrack*>(ms);
      }
      if (rot) {
        const Identifier& id = rot->identify();
        if (m_doSharedHits and m_usePixel and m_pixelId->is_pixel(id)) {
          // check if shared
          bool hitIsSplit(false);
          if (m_runningTIDE_Ambi) {
            const PixelClusterOnTrack* pix = nullptr;
            if (rot->rioType(Trk::RIO_OnTrackType::PixelCluster)) {
              pix = static_cast<const PixelClusterOnTrack*>(rot);
            }
            if (pix) {
              const InDet::PixelCluster* pixPrd = pix->prepRawData();
              const Trk::ClusterSplitProbabilityContainer::ProbabilityInfo&
                splitProb = getClusterSplittingProbability(ctx, pixPrd);
              if (pixPrd and splitProb.isSplit()) {
                summary.m_information[Trk::numberOfPixelSplitHits]++;
                hitIsSplit = true;
                if (m_pixelId->is_barrel(id) and
                    m_pixelId->layer_disk(id) == 0) {
                  summary.m_information[Trk::numberOfInnermostLayerSplitHits]++;
                }
                if (m_pixelId->is_barrel(id) and
                    m_pixelId->layer_disk(id) == 1) {
                  summary.m_information
                    [Trk::numberOfNextToInnermostLayerSplitHits]++;
                }
              }
            }
          }
          // If we are running the TIDE ambi don't count split hits as shared
          if (not(m_runningTIDE_Ambi and hitIsSplit)) {
            if (isShared(prd_to_track_map, m_assoTool, *(rot->prepRawData()))) {
              ATH_MSG_DEBUG("shared Pixel hit found");
              summary.m_information[Trk::numberOfPixelSharedHits]++;
              if ((m_pixelId->is_barrel(id) and
                   m_pixelId->layer_disk(id) == 0)) {
                ATH_MSG_DEBUG(
                  "--> shared Pixel hit is in Innermost Pixel layer");
                summary
                  .m_information[Trk::numberOfInnermostPixelLayerSharedHits]++;
              } else if ((m_pixelId->is_barrel(id) and
                          m_pixelId->layer_disk(id) == 1)) {
                ATH_MSG_DEBUG(
                  "--> shared Pixel hit is in Next To Innermost Pixel layer");
                summary.m_information
                  [Trk::numberOfNextToInnermostPixelLayerSharedHits]++;
              }
            }
          }
        } else if (m_doSharedHits and m_useSCT and m_sctId->is_sct(id)) {
          // used in more than one track ?
          if (isShared(prd_to_track_map, m_assoTool, *(rot->prepRawData()))) {
            ATH_MSG_DEBUG("shared SCT hit found");
            summary.m_information[Trk::numberOfSCTSharedHits]++;
          }
        }
        if (m_doSharedHitsTRT and m_useTRT and m_trtId->is_trt(id)) {
          // used in more than one track ?
          if (isShared(prd_to_track_map, m_assoTool, *(rot->prepRawData()))) {
            ATH_MSG_DEBUG("shared TRT hit found");
            summary.m_information[Trk::numberOfTRTSharedHits]++;
          }
        }
      }
    }
  }
  }

void
InDet::InDetTrackSummaryHelperTool::addDetailedTrackSummary(
  const EventContext&,
  const Trk::Track&,
  Trk::TrackSummary&) const
{
}

StatusCode
InDet::InDetTrackSummaryHelperTool::finalize()
{
  return StatusCode::SUCCESS;
}
