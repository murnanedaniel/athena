/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/***************************************************************************
                          Summary.cxx  -  description
                             -------------------
    begin                : Tue Feb 17 2004
    email                : edward.moyse@cern.ch
 ***************************************************************************/

#include "TrkTrackSummary/TrackSummary.h"
#include "GaudiKernel/MsgStream.h"

std::atomic<unsigned int> Trk::TrackSummary::s_numberOfInstantiations{ 0 };
const int Trk::TrackSummary::SummaryTypeNotSet = -1;

Trk::TrackSummary::TrackSummary()
  : m_information(numberOfTrackSummaryTypes, SummaryTypeNotSet)
  , m_idHitPattern(0)
  , m_muonTrackSummary(nullptr)
{
#ifndef NDEBUG
  s_numberOfInstantiations++; // new TrackSummary, so increment total count
#endif
}

Trk::TrackSummary::TrackSummary(const std::vector<int>& information,
                                std::bitset<numberOfDetectorTypes>& hitPattern)
  : m_information(information)
  , m_idHitPattern(hitPattern.to_ulong())
  , m_muonTrackSummary(nullptr)
{
#ifndef NDEBUG
  s_numberOfInstantiations++; // new TrackSummary, so increment total count
#endif
}

Trk::TrackSummary::TrackSummary(const TrackSummary& rhs)
  : m_information(rhs.m_information)
  , m_idHitPattern(rhs.m_idHitPattern)
{
#ifndef NDEBUG
  s_numberOfInstantiations++; // new TrackSummary, so increment total count
#endif
  if (rhs.m_muonTrackSummary) {
    m_muonTrackSummary =
      std::make_unique<MuonTrackSummary>(*rhs.m_muonTrackSummary);
  } else {
    m_muonTrackSummary = nullptr;
  }
}

Trk::TrackSummary&
Trk::TrackSummary::operator=(const TrackSummary& rhs)
{
  if (&rhs != this) {
    m_information = rhs.m_information;
    m_idHitPattern = rhs.m_idHitPattern;
    m_muonTrackSummary.reset(rhs.m_muonTrackSummary
                               ? new MuonTrackSummary(*rhs.m_muonTrackSummary)
                               : nullptr);
  }
  return *this;
}

Trk::TrackSummary::~TrackSummary()
{
#ifndef NDEBUG
  s_numberOfInstantiations--; // delete TrackSummary, so decrement total count
#endif
}

Trk::TrackSummary&
Trk::TrackSummary::operator+=(const TrackSummary& ts)
{
  if (this != &ts) {
    for (int i = 0; i < numberOfTrackSummaryTypes; ++i) {
      // if added number is <0, leave as is (this also catches the case where
      // both are <0)
      if (ts.m_information[i] < 0) {
        continue;
      }

      if (m_information[i] < 0) {
        m_information[i]++;
      }
      m_information[i] += ts.m_information[i];
    }
    if (!m_muonTrackSummary) {
      m_muonTrackSummary.reset(ts.m_muonTrackSummary
                                 ? new MuonTrackSummary(*ts.m_muonTrackSummary)
                                 : nullptr);
    }
  }
  return *this;
}

namespace Trk {
namespace {

template<class T_out>
inline T_out&
dumpTrackSummary(T_out& out, const TrackSummary& trackSum)
{
  out << "Persistant track summary information:"
      << "\n";
  out << " * Number of contrib. Pixel Layer: "
      << trackSum.get(numberOfContribPixelLayers) << "\n";
  out << " * Number of Innermost Pixel layer hits        : "
      << trackSum.get(numberOfInnermostPixelLayerHits) << "\n";
  out << " * Number of Innermost Pixel layer shared hits : "
      << trackSum.get(numberOfInnermostPixelLayerSharedHits) << "\n";
  out << " * Number of Innermost Pixel layer outliers    : "
      << trackSum.get(numberOfInnermostPixelLayerOutliers) << "\n";
  out << " * Number of Next-To-Innermost Pixel layer hits        : "
      << trackSum.get(numberOfNextToInnermostPixelLayerHits) << "\n";
  out << " * Number of Next-To-Innermost Pixel layer shared hits : "
      << trackSum.get(numberOfNextToInnermostPixelLayerSharedHits) << "\n";
  out << " * Number of Next-To-Innermost Pixel layer outliers    : "
      << trackSum.get(numberOfNextToInnermostPixelLayerOutliers) << "\n";
  out << " * Number of pixel hits          : "
      << trackSum.get(numberOfPixelHits) << "\n";
  out << " * Number of pixel outliers      : "
      << trackSum.get(numberOfPixelOutliers) << "\n";
  out << " * Number of spoilt pixel hits   : "
      << trackSum.get(numberOfPixelSpoiltHits) << "\n";
  out << " * Number of pixel holes         : "
      << trackSum.get(numberOfPixelHoles) << "\n";
  out << " * Number of pixel shared hits   : "
      << trackSum.get(numberOfPixelSharedHits) << "\n";
  out << " * Number of GangedPixels        : "
      << trackSum.get(numberOfGangedPixels) << "\n";
  out << " * Number of GangedFlaggedFakes  : "
      << trackSum.get(numberOfGangedFlaggedFakes) << "\n";
  out << " * Number of dead pixel sensors  : "
      << trackSum.get(numberOfPixelDeadSensors) << "\n";
  out << " * Number of DBM hits            : " << trackSum.get(numberOfDBMHits)
      << "\n";
  out << " * Number of SCT hits            : " << trackSum.get(numberOfSCTHits)
      << "\n";
  out << " * Number of SCT outliers        : "
      << trackSum.get(numberOfSCTOutliers) << "\n";
  out << " * Number of SCT holes           : " << trackSum.get(numberOfSCTHoles)
      << "\n";
  out << " * Number of SCT double Holes    : "
      << trackSum.get(numberOfSCTDoubleHoles) << "\n";
  out << " * Number of SCT shared hits     : "
      << trackSum.get(numberOfSCTSharedHits) << "\n";
  out << " * Number of dead SCT sensors    : "
      << trackSum.get(numberOfSCTDeadSensors) << "\n";
  out << " * Number of spoilt SCT hits     : "
      << trackSum.get(numberOfSCTSpoiltHits) << "\n";
  out << " * Number of TRT hits            : " << trackSum.get(numberOfTRTHits)
      << "\n";
  out << " * Number of TRT xenon hits      : "
      << trackSum.get(numberOfTRTXenonHits) << "\n";
  out << " * Number of TRT outliers        : "
      << trackSum.get(numberOfTRTOutliers) << "\n";
  out << " * Number of TRT TR hits         : "
      << trackSum.get(numberOfTRTHighThresholdHits) << "\n";
  out << " * Number of TRT TR hits total   : "
      << trackSum.get(numberOfTRTHighThresholdHitsTotal) << "\n";
  out << " * Number of TRT TR outliers     : "
      << trackSum.get(numberOfTRTHighThresholdOutliers) << "\n";
  out << " * Number of TRT holes           : " << trackSum.get(numberOfTRTHoles)
      << "\n";
  out << " * Number of TRT tube hits       : "
      << trackSum.get(numberOfTRTTubeHits) << "\n";
  out << " * Number of TRT Shared hits     : "
      << trackSum.get(numberOfTRTSharedHits) << "\n";
  out << " * Number of dead TRT straws     : "
      << trackSum.get(numberOfTRTDeadStraws) << "\n";
  out << " * Number of MDT hits            : " << trackSum.get(numberOfMdtHits)
      << "\n";
  out << " * Number of MDT holes           : " << trackSum.get(numberOfMdtHoles)
      << "\n";

  // New Small Wheel
  out << " * Number of STGC hits           : phi "
      << trackSum.get(numberOfStgcPhiHits) << " eta "
      << trackSum.get(numberOfStgcEtaHits) << "\n";
  out << " * Number of STGC holes          : phi "
      << trackSum.get(numberOfStgcPhiHoles) << " eta  "
      << trackSum.get(numberOfStgcEtaHoles) << "\n";
  out << " * Number of MM hits             : " << trackSum.get(numberOfMmHits)
      << "\n";
  out << " * Number of MM holes            : " << trackSum.get(numberOfMmHoles)
      << "\n";

  out << " * Number of TGC hits            : phi "
      << trackSum.get(numberOfTgcPhiHits) << " eta "
      << trackSum.get(numberOfTgcEtaHits) << "\n";
  out << " * Number of TGC holes           : phi "
      << trackSum.get(numberOfTgcPhiHoles) << " eta  "
      << trackSum.get(numberOfTgcEtaHoles) << "\n";
  out << " * Number of RPC hits            : phi "
      << trackSum.get(numberOfRpcPhiHits) << " eta "
      << trackSum.get(numberOfRpcEtaHits) << "\n";
  out << " * Number of RPC holes           : phi "
      << trackSum.get(numberOfRpcPhiHoles) << " eta  "
      << trackSum.get(numberOfRpcEtaHoles) << "\n";
  out << " * Number of CSC hits            : phi "
      << trackSum.get(numberOfCscPhiHits) << " eta "
      << trackSum.get(numberOfCscEtaHits) << "\n";
  out << " * Number of unspoilt CSC etahits:"
      << trackSum.get(numberOfCscUnspoiltEtaHits) << "\n";
  out << " * Number of CSC holes           : phi "
      << trackSum.get(numberOfCscPhiHoles) << " eta  "
      << trackSum.get(numberOfCscEtaHoles) << "\n";
  out << " * Number of Outliers            : "
      << trackSum.get(numberOfOutliersOnTrack) << "\n";
  out << "Transient track summary information:"
      << "\n";
  out << " * scatter of chi2 on Surface    : "
      << float(trackSum.get(standardDeviationOfChi2OS)) / 100. << "\n";
  // this is a bit nasty, but I don't have access to internal data members
  out << " Hit pattern (see DetectorType enum for meaning) : ";
  for (int i = 0; i < Trk::numberOfDetectorTypes; ++i) {
    if (trackSum.isHit(static_cast<Trk::DetectorType>(i))) {
      out << "1";
    } else {
      out << "0";
    }
  }
  out << "\n";
  if (nullptr != trackSum.muonTrackSummary()) {
    out << *(trackSum.muonTrackSummary());
  }
  return out;
}
}
}

std::ostream&
Trk::operator<<(std::ostream& out, const TrackSummary& trackSum)
{
  return dumpTrackSummary(out, trackSum);
}

MsgStream&
Trk::operator<<(MsgStream& out, const TrackSummary& trackSum)
{
  return dumpTrackSummary(out, trackSum);
}

unsigned int
Trk::TrackSummary::numberOfInstantiations()
{
  return s_numberOfInstantiations;
}

