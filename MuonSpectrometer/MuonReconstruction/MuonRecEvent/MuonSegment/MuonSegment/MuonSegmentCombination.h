/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// MdtSegment.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef MUONSEGMENT_MUONSEGMENTCOMBINATION_H
#define MUONSEGMENT_MUONSEGMENTCOMBINATION_H

namespace Trk {
  class Segment;
}

namespace Muon {
  class MuonSegment;
}

#include <vector>


namespace Muon {

  /** @class MuonSegmentCombination
      
  Class to hold a set of MuonSegments belonging together.
  The MuonSegments are stored per chamber, making it easier to keep track of possible
  ambiguities per station.
  
  @author niels.van.eldik@cern.ch
  */
  class MuonSegmentCombination {
  public:
    typedef std::vector< const MuonSegment* > SegmentVec;
    typedef std::vector< SegmentVec* >    SegmentVecVec;
  public:
    /** Default constructor */
    MuonSegmentCombination();
    
    /** Destructor */
    ~MuonSegmentCombination();

    /** Copy constructor */
    MuonSegmentCombination( const MuonSegmentCombination& );

    /** assigment operator */
    MuonSegmentCombination& operator=( const MuonSegmentCombination& );

    /** Add a set of Segments for a give station.
	For now no sorting but this could be changed so the segments
	are stored with increasing radius. Also no check is performed 
	whether there are already segments for the given station. 
	This is up to the user.
    */
    bool addSegments( SegmentVec* );

    /** Number of stations with segment */
    unsigned int numberOfStations() const;

    /** Access to segments in a given station */
    const SegmentVec* stationSegments( unsigned int index ) const;

    /** Number of ambiguities */
    unsigned int numberOfAmbiguities() const;

    void setUse2LayerSegments(bool use2Lay){use2LayerSegs=use2Lay;}

    bool use2LayerSegments() const;

  private:
    /** clear data */
    void clear();

    /** copy data */
    void copy( const MuonSegmentCombination& segc );

    SegmentVecVec         m_segmentsPerStation;

    //if the station is a CSC station with 2-layer segment finding enabled
    bool use2LayerSegs;
  };

  inline bool MuonSegmentCombination::use2LayerSegments() const
    {
      return use2LayerSegs;
    }

  inline  bool MuonSegmentCombination::addSegments( SegmentVec* segs )
  {
    m_segmentsPerStation.push_back( segs );
    return true;
  }

  inline unsigned int MuonSegmentCombination::numberOfStations() const 
  {
    return m_segmentsPerStation.size();
  }
  
  inline const MuonSegmentCombination::SegmentVec* 
    MuonSegmentCombination::stationSegments( unsigned int index ) const
  {
    if( index < numberOfStations() ) return m_segmentsPerStation[index];
    return 0;
  }

  inline unsigned int MuonSegmentCombination::numberOfAmbiguities() const
  {
    unsigned int solutions(1);
    SegmentVecVec::const_iterator it = m_segmentsPerStation.begin();
    SegmentVecVec::const_iterator it_end = m_segmentsPerStation.end();
    for( ; it!=it_end; ++it ) solutions *= (*it)->size();
    return solutions;
  }

}

#endif
