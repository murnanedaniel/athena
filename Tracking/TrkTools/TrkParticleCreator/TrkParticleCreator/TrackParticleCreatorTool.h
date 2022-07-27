/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************
TrackParticleCreatorTool.h  -  Description
-------------------
begin   : Autumn 2003
authors : Andreas Wildauer (CERN PH-ATC), Fredrik Akesson (CERN PH-ATC)
email   : andreas.wildauer@cern.ch, fredrik.akesson@cern.ch
changes : 11.02.04 added docu

***************************************************************************/
#ifndef TRKPARTICLECREATOR_PARTICLECREATORTOOL_H
#define TRKPARTICLECREATOR_PARTICLECREATORTOOL_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "TrkToolInterfaces/ITrackParticleCreatorTool.h"

#include "AthContainers/AuxElement.h"
#include "EventPrimitives/EventPrimitivesHelpers.h"
#include "ITrackToVertex/ITrackToVertex.h"
#include "InDetIdentifier/PixelID.h"
#include "MuonRecToolInterfaces/IMuonHitSummaryTool.h"
#include "PixelGeoModel/IBLParameterSvc.h"
#include "TrkEventPrimitives/FitQuality.h"
#include "TrkExInterfaces/IExtrapolator.h"
#include "TrkParameters/TrackParameters.h"
#include "TrkParticleBase/TrackParticleBase.h" // for TrackParticleOrigin enum
#include "TrkToolInterfaces/IExtendedTrackSummaryTool.h"
#include "TrkToolInterfaces/IPixelToTPIDTool.h" //template parameter to tool handle
#include "TrkToolInterfaces/ITRT_ElectronPidTool.h" //template parameter to tool handle
#include "InDetRecToolInterfaces/IInDetTestPixelLayerTool.h"
#include "TrkTrack/TrackInfo.h"
#include "TrkTrackSummary/MuonTrackSummary.h"
#include "TrkTrackSummary/TrackSummary.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackingPrimitives.h"
#include "xAODTracking/VertexFwd.h"

#include "GaudiKernel/ToolHandle.h"
// MagField cache
#include "MagFieldConditions/AtlasFieldCacheCondObj.h"
#include "MagFieldElements/AtlasFieldCache.h"

namespace Rec {
class TrackParticle;
}

namespace InDet {
class BeamSpotData;
}

namespace Trk {

class Track;
class VxCandidate;

class TrackParticleCreatorTool final
  : public extends<AthAlgTool, ITrackParticleCreatorTool>
{
public:
  TrackParticleCreatorTool(const std::string&,
                           const std::string&,
                           const IInterface*);

  virtual ~TrackParticleCreatorTool() = default;

  virtual StatusCode initialize() override;

  /** Method to construct a xAOD::TrackParticle from a Rec::TrackParticle.
  @param track particle
  @param TrackParticleContainer needed to have an AuxStore, if provided particle
  will be added to store which takes ownership
  */
  virtual xAOD::TrackParticle* createParticle(
    const EventContext& ctx,
    const Rec::TrackParticle& trackParticle,
    xAOD::TrackParticleContainer* container) const override final;

  /** Method to construct a xAOD::TrackParticle from a passed Track.
  Will keep parameters  based on m_keepParameters,m_keepFirstParameters,
  m_keepAllPerigee.
  It will use the exising summary or redo it based on m_useTrackSummaryTool
  @param track Pointer to a valid track (i.e. do not pass a zero!). Ownership is
  not taken (i.e. it will not be deleted)
  @param TrackParticleContainer needed to have an AuxStore, if provided particle
  will be added to store which takes ownership
  @param xAOD::Vertex Pointer to a valid vxCandidate (i.e. do not pass a zero!).
  Ownership is not taken (i.e. it will not be deleted)
  @param prtOrigin Particle type
  @param prd_to_track_map an optional PRD-to-track map to compute shared hits.
  */
  virtual xAOD::TrackParticle* createParticle(
    const EventContext& ctx,
    const Trk::Track& track,
    xAOD::TrackParticleContainer* container,
    const xAOD::Vertex* vxCandidate,
    xAOD::ParticleHypothesis prtOrigin,
    const Trk::PRDtoTrackMap* prd_to_track_map) const override final;

  /** Method to construct a TrackParticle from a passed Track.
  Will keep parameters  based on m_keepParameters,m_keepFirstParameters,
  m_keepAllPerigee.
  It will use the exising summary or redo it based on m_useTrackSummaryTool
  @param track element link to a valid track (i.e. do not pass a zero!).
  @param TrackParticleContainer needed to have an AuxStore, if provided particle
  will be added to store which takes ownership
  @param xAOD::Vertex Pointer to a valid vxCandidate (i.e. do not pass a zero!).
  Ownership is not taken (i.e. it will not be deleted)
  @param prtOrigin
  @param prd_to_track_map an optional PRD-to-track map to compute shared hits.
  */
  virtual xAOD::TrackParticle* createParticle(
    const EventContext& ctx,
    const ElementLink<TrackCollection>& trackLink,
    xAOD::TrackParticleContainer* container,
    const xAOD::Vertex* vxCandidate,
    xAOD::ParticleHypothesis prtOrigin,
    const Trk::PRDtoTrackMap* prd_to_track_map) const override final;

  /** create a xAOD::TrackParticle out of constituents */
  virtual xAOD::TrackParticle* createParticle(
    const EventContext& ctx,
    const Perigee* perigee,
    const FitQuality* fq,
    const TrackInfo* trackInfo,
    const TrackSummary* summary,
    const std::vector<const Trk::TrackParameters*>& parameters,
    const std::vector<xAOD::ParameterPosition>& positions,
    xAOD::ParticleHypothesis prtOrigin,
    xAOD::TrackParticleContainer* container,
    bool addInfoIfMuon) const override final;

  virtual const InDet::BeamSpotData* CacheBeamSpotData(
    const EventContext& ctx) const override final;

  /** Method to set FitQuality of a xAOD::TrackParticle */
  void setFitQuality(xAOD::TrackParticle& tp, const FitQuality& fq) const;

  /** Method to set TrackInfo of a xAOD::TrackParticle */
  void setTrackInfo(xAOD::TrackParticle& tp,
                    const TrackInfo& trackInfo,
                    xAOD::ParticleHypothesis prtOrigin) const;

  /** Method to set TrackSummary of a xAOD::TrackParticle */
  void setTrackSummary(xAOD::TrackParticle& tp,
                       const TrackSummary& summary) const;

  /** Add Pixel and TRT PID information to the track particle.
   * @param ctx the current event context
   * @param track a valid track or nullptr to set all PID values to the default value.
   * @param tp the trackparticle in which the PID values are set.
   */
   void addPIDInformation(const EventContext& ctx, const Track *track, xAOD::TrackParticle& tp) const;

   /** Add extra detailed hit summary info not computed in Trk::TrkSummary */
   void addDetailedHitInformation(const DataVector<const TrackStateOnSurface>* trackStates, xAOD::TrackParticle& tp) const;

   /** Add expected hit info for innermost pixel layers not computed in Trk::TrkSummary */
   void addExpectedHitInformation(const Perigee* perigee, xAOD::TrackParticle& tp) const;

  /** Method to set Defining parameters of a xAOD::TrackParticle */
  void setDefiningParameters(xAOD::TrackParticle& tp,
                             const Perigee& perigee) const;

  /** Method to set parameters of a xAOD::TrackParticle */
  void setParameters(
    const EventContext& ctx,
    xAOD::TrackParticle& tp,
    const std::vector<const Trk::TrackParameters*>& parameters,
    const std::vector<xAOD::ParameterPosition>& positions) const;

  static void setTilt(xAOD::TrackParticle& tp, float tiltx, float tilty);

  static void setHitPattern(xAOD::TrackParticle& tp, unsigned long hitpattern);

  static void setNumberOfUsedHits(xAOD::TrackParticle& tp, int hits);

  static void setNumberOfOverflowHits(xAOD::TrackParticle& tp, int overflows);

  /** Get the name used for the decoration of the track particle with the number
   * of used hits for TRT dE/dx computation.*/
  static const std::string& trtdEdxUsedHitsAuxName()
  {
    return s_trtdEdxUsedHitsDecorationName;
  }

protected:
  /** create a xAOD::TrackParticle out of constituents */
  xAOD::TrackParticle* createParticle(
    const EventContext& ctx,
    const Perigee* perigee,
    const FitQuality* fq,
    const TrackInfo* trackInfo,
    const TrackSummary* summary,
    const std::vector<const Trk::TrackParameters*>& parameters,
    const std::vector<xAOD::ParameterPosition>& positions,
    xAOD::ParticleHypothesis prtOrigin,
    xAOD::TrackParticleContainer* container,
    const Trk::Track *track,
    bool addInfoIfMuon = false) const;

private:
  void compare(const Rec::TrackParticle& tp,
               const xAOD::TrackParticle& tpx) const;
  void compare(const TrackParameters& tp1, const TrackParameters& tp2) const;
  /**atlas id helper*/
  const AtlasDetectorID* m_detID;
  const PixelID* m_pixelID;

  // Need to change to private when is safe to do so
  PublicToolHandle<IExtendedTrackSummaryTool> m_trackSummaryTool{
    this,
    "TrackSummaryTool",
    "Trk::TrackSummaryTool/AtlasTrackSummaryTool"
  };

  ToolHandle<Reco::ITrackToVertex> m_trackToVertex{
    this,
    "TrackToVertex",
    "Reco::TrackToVertex/TrackToVertex"
  };
  ToolHandle<Muon::IMuonHitSummaryTool> m_hitSummaryTool{
    this,
    "MuonSummaryTool",
    ""
  };

  ServiceHandle<IBLParameterSvc> m_IBLParameterSvc{
    this,
    "IBLParameterSvc",
    "IBLParameterSvc"
  };

  SG::ReadCondHandleKey<AtlasFieldCacheCondObj> m_fieldCacheCondObjInputKey{
    this,
    "AtlasFieldCacheCondObj",
    "fieldCondObj",
    "Name of the Magnetic Field conditions object key"
  };

   /**tool to calculate electron probabilities*/
  ToolHandle<ITRT_ElectronPidTool> m_eProbabilityTool{ this,
                                                       "TRT_ElectronPidTool",
                                                       "",
                                                       "" };
  /**tool to calculate dE/dx using pixel clusters*/
  ToolHandle<IPixelToTPIDTool> m_dedxtool{ this, "PixelToTPIDTool", "", "" };

  /**tool to calculate expected hit information in innermost layers*/
  ToolHandle<InDet::IInDetTestPixelLayerTool> m_testPixelLayerTool{ this,
                                                                    "TestPixelLayerTool",
                                                                    "",
                                                                    "" };

  /** Configurable to set the eProbabilities and extra track summary types which
   * are to be copied from the track summary.*/
  std::vector<std::string> m_copyExtraSummaryName;

  /** Enums of an eProbability which are set in the xAOD::TrackSummary.*/
  std::vector<Trk::eProbabilityType> m_copyEProbabilities;

  /** The pairs if enums  of an eProbability which is added as a decoration to
   * the track particle and the name of the decoration.*/
  std::vector<std::pair<SG::AuxElement::Accessor<float>, Trk::eProbabilityType>>
    m_decorateEProbabilities;
  std::vector<std::pair<SG::AuxElement::Accessor<uint8_t>, Trk::SummaryType>>
    m_decorateSummaryTypes;

  /** Name used for the decoration of the track particle with TRT dE/dx .*/
  static const std::string s_trtdEdxUsedHitsDecorationName;
  static const SG::AuxElement::Accessor<uint8_t> s_trtdEdxUsedHitsDecoration;

  bool m_doIBL;
  bool m_doITk;
  ///< if the track contains a summary, the shared, expected hit, and PID
  ///< information will be recomputed. The summary of the track is not updated.
  bool m_computeAdditionalInfo;
  /** the following keep options are mutually exclusive **/
  /** keep all TrackParameters */
  bool m_keepParameters;
  ///< keep the first parameters when creating track particles.
  bool m_keepFirstParameters;
  /** keep all MeasuredPerigee parameters (e.g. adding those that may exist at
   * Volume boundaries) */
  bool m_keepAllPerigee;
  bool m_expressPerigeeToBeamSpot;
  int m_badclusterID;

  std::string m_perigeeExpression;
  std::vector<std::string> m_perigeeOptions{ "BeamLine",
                                             "BeamSpot",
                                             "Vertex",
                                             "Origin" };

  bool m_checkConversion;
  int m_minSiHits;
  double m_minPt;
};

} // end of namespace Trk
#include "TrkParticleCreator/TrackParticleCreatorTool.icc"
#endif
