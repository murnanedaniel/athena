/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


#ifndef INDETCOSMICSCORINGTOOL_H
#define INDETCOSMICSCORINGTOOL_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrkEventPrimitives/TrackScore.h"
#include "TrkToolInterfaces/ITrackScoringTool.h"
#include "TrkToolInterfaces/ITrackSummaryTool.h"
#include <vector>
#include <string>

namespace Trk {
  class IExtrapolator;
  class Track;
  class TrackSummary;
}


namespace InDet {


/**Concrete implementation of the ITrackScoringTool pABC*/
class InDetCosmicScoringTool : virtual public Trk::ITrackScoringTool, public AthAlgTool
{

public:
  InDetCosmicScoringTool(const std::string&,const std::string&,const IInterface*);
  virtual ~InDetCosmicScoringTool ();
  virtual StatusCode initialize();
  virtual StatusCode finalize  ();
  /** create a score based on how good the passed track is*/
  Trk::TrackScore score( const Trk::Track& track, const bool suppressHoleSearch );
  
  /** create a score based on how good the passed TrackSummary is*/
  Trk::TrackScore simpleScore( const Trk::Track& track, const Trk::TrackSummary& trackSum );
  
 private:
  
  /**\todo make this const, once createSummary method is const*/
  ToolHandle<Trk::ITrackSummaryTool> m_trkSummaryTool;

  int m_nWeightedClustersMin; 
  int m_minTRTHits;

};


}
#endif 
