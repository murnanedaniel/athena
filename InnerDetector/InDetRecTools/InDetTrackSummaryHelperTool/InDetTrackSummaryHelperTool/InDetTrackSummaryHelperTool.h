/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef INDETTRACKSUMMARYHELPERTOOL_H
#define INDETTRACKSUMMARYHELPERTOOL_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "TrkToolInterfaces/IExtendedTrackSummaryHelperTool.h"

#include "TRT_ConditionsServices/ITRT_StrawStatusSummaryTool.h"
#include "TrkEventPrimitives/ParticleHypothesis.h"
#include "TrkToolInterfaces/IPRD_AssociationTool.h"
#include "TrkToolInterfaces/IPixelToTPIDTool.h"
#include "TrkToolInterfaces/ITrackHoleSearchTool.h"
#include "TrkTrackSummary/TrackSummary.h" // defines the Trk::numberOfDetectorTypes enum

#include "InDetPrepRawData/PixelCluster.h"
#include "TrkEventUtils/ClusterSplitProbabilityContainer.h"

#include "GaudiKernel/ToolHandle.h"

#include <bitset>
#include <vector>

class ITRT_StrawSummaryTool;
class PixelID;
class SCT_ID;
class TRT_ID;

namespace Trk {
class CompetingRIOsOnTrack;
class RIO_OnTrack;
class Track;
class TrackStateOnSurface;
}

namespace InDet {

class InDetTrackSummaryHelperTool final
  : public extends<AthAlgTool, Trk::IExtendedTrackSummaryHelperTool>
{
public:
  /** constructor */
  InDetTrackSummaryHelperTool(const std::string&,
                              const std::string&,
                              const IInterface*);

  /** destructor */
  virtual ~InDetTrackSummaryHelperTool() = default;

  /** standard AlgTool methods: initialise retrieves Tools, finalize does
   * nothing */
  virtual StatusCode initialize() override;
  virtual StatusCode finalize() override;

  /** Input : rot, tsos
      Output: Changes in information and hitPattern
      Input quantities rot, tsos are used to increment the counts for hits and
     outliers in information and to set the proper bits in hitPattern.
  */
  using IExtendedTrackSummaryHelperTool::addDetailedTrackSummary;
  using IExtendedTrackSummaryHelperTool::analyse;

  virtual void analyse(
    const EventContext& ctx,
    const Trk::Track& track,
    const Trk::PRDtoTrackMap* prd_to_track_map,
    const Trk::RIO_OnTrack* rot,
    const Trk::TrackStateOnSurface* tsos,
    std::vector<int>& information,
    std::bitset<Trk::numberOfDetectorTypes>& hitPattern) const override final;

  virtual void analyse(
    const EventContext& ctx,
    const Trk::Track& track,
    const Trk::PRDtoTrackMap* prd_to_track_map,
    const Trk::CompetingRIOsOnTrack* crot,
    const Trk::TrackStateOnSurface* tsos,
    std::vector<int>& information,
    std::bitset<Trk::numberOfDetectorTypes>& hitPattern) const override final;
  /** @copydoc Trk::ITrackSummaryHelperTool::addDetailedTrackSummary(const
   * Trk::Track&, Trk::TrackSummary&)*/

  virtual void addDetailedTrackSummary(const EventContext& ctx,
                                       const Trk::Track&,
                                       Trk::TrackSummary&) const override final;

  /** Input : track, partHyp
      Output: Changes in information
      This method first calls the method getListOfHits to isolate the relevant
     hits on the track before calling the method performHoleSearchStepWise
     which then performs the actual hole search. Additionally the Layers of
     the Pixel Detector which contribute measurements to the track are counted
      If problems occur, the information counters for Holes and PixelLayers
     are reset to -1 flagging them as not set.
  */
  virtual void searchForHoles(
    const Trk::Track& track,
    std::vector<int>& information,
    const Trk::ParticleHypothesis partHyp = Trk::pion) const override final;

  /** this method simply updaes the shared hit content - it is
   * designed/optimised for track collection merging */
  virtual void updateSharedHitCount(
    const Trk::Track& track,
    const Trk::PRDtoTrackMap* prd_to_track_map,
    Trk::TrackSummary& summary) const override final;

private:
  const Trk::ClusterSplitProbabilityContainer::ProbabilityInfo&
  getClusterSplittingProbability(const EventContext& ctx,
                                 const InDet::PixelCluster* pix) const
  {
    if (!pix || m_clusterSplitProbContainer.key().empty()) {
      return Trk::ClusterSplitProbabilityContainer::getNoSplitProbability();
    }
    SG::ReadHandle<Trk::ClusterSplitProbabilityContainer> splitProbContainer(
      m_clusterSplitProbContainer, ctx);
    if (!splitProbContainer.isValid()) {
      ATH_MSG_FATAL("Failed to get cluster splitting probability container "
                    << m_clusterSplitProbContainer);
    }
    return splitProbContainer->splitProbability(pix);
  }
  /**ID pixel helper*/
  const PixelID* m_pixelId{ nullptr };

  /**ID SCT helper*/
  const SCT_ID* m_sctId{ nullptr };

  /**ID TRT helper*/
  const TRT_ID* m_trtId{ nullptr };

  /**Association tool - used to work out which (if any)
   * PRDs are shared between tracks*/
  PublicToolHandle<Trk::IPRD_AssociationTool> m_assoTool{
    this,
    "AssoTool",
    "InDet::InDetPRD_AssociationToolGangedPixels"
  };
  PublicToolHandle<Trk::ITrackHoleSearchTool> m_holeSearchTool{
    this,
    "HoleSearch",
    "InDet::InDetTrackHoleSearchTool"
  };
  ToolHandle<Trk::IPixelToTPIDTool> m_pixeldedxtool{ this,
                                                     "PixelToTPIDTool",
                                                     "" };
  ToolHandle<ITRT_StrawStatusSummaryTool> m_TRTStrawSummaryTool{
    this,
    "TRTStrawSummarySvc",
    "TRT_StrawStatusSummaryTool",
    "The ConditionsSummaryTool"
  };

  SG::ReadHandleKey<Trk::ClusterSplitProbabilityContainer>
    m_clusterSplitProbContainer{ this, "ClusterSplitProbabilityName", "", "" };

  BooleanProperty m_usePixel{ this, "usePixel", true };
  BooleanProperty m_useSCT{ this, "useSCT", true };
  BooleanProperty m_useTRT{ this, "useTRT", true };
  BooleanProperty m_doSharedHits{ this, "DoSharedHits", false };
  BooleanProperty m_doSharedHitsTRT{ this, "DoSharedHitsTRT", false };
  BooleanProperty m_doSplitPixelHits{ this, "DoSplitHits", true };
  BooleanProperty m_overwriteidsummary{ this, "OverwriteIDSummary", false };
  BooleanProperty m_runningTIDE_Ambi{ this, "RunningTIDE_Ambi", false };
};

}
#endif
