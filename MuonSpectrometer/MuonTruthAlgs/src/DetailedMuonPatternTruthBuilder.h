/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DETAILEDMUONPATTERNTRUTHBUILDER_H
#define DETAILEDMUONPATTERNTRUTHBUILDER_H

#include "MuonRecToolInterfaces/IDetailedMuonPatternTruthBuilder.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

#include "AthLinks/ElementLink.h"
#include "AtlasDetDescr/AtlasDetectorID.h"
#include "TrkTruthData/PRD_MultiTruthCollection.h"
#include "TrkEventUtils/InverseMultiMap.h"
#include "MuonPattern/MuonPatternCollection.h"
#include "MuonPattern/MuonPatternCombination.h"
#include "MuonPattern/MuonPatternCombinationCollection.h"
#include "MuonPattern/DetailedMuonPatternTruthCollection.h"
#include "TrkToolInterfaces/ITruthTrajectoryBuilder.h"
#include "TrkTruthData/DetailedSegmentTruth.h"
#include "TrkTruthData/SubDetHitStatistics.h"
#include "TrkTruthData/DetailedTrackTruthCollection.h"
#include "TrkTruthData/DetailedTrackTruth.h"
#include "TrkTruthData/TruthTrajectory.h"
#include "TrkTrack/Track.h"
#include "TrkRIO_OnTrack/RIO_OnTrack.h"
#include "TrkCompetingRIOsOnTrack/CompetingRIOsOnTrack.h"
#include "TrkMeasurementBase/MeasurementBase.h"
#include "TrkPrepRawData/PrepRawData.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonSimData/MuonSimDataCollection.h"
#include "MuonSimData/CscSimDataCollection.h"
#include "TrackRecord/TrackRecord.h"
#include "TrkSegment/SegmentCollection.h"
#include "TrkParameters/TrackParameters.h"
#include "MuonRecToolInterfaces/IMdtDriftCircleOnTrackCreator.h"
#include "MuonRecToolInterfaces/IMuonClusterOnTrackCreator.h"

#include <string>
#include <vector>

namespace Trk {
  
  class DetailedMuonPatternTruthBuilder: virtual public IDetailedMuonPatternTruthBuilder, public AthAlgTool
  {
  public:
    DetailedMuonPatternTruthBuilder(const std::string& type, const std::string& name, const IInterface* parent);

    virtual StatusCode initialize();

    /** See description for IDetailedMuonPatternTruthBuilder::buildDetailedTrackTruth() */
    virtual void buildDetailedMuonPatternTruth(DetailedMuonPatternTruthCollection *output,
					 const MuonPatternCombinationCollection& tracks,
					 const std::vector<const PRD_MultiTruthCollection*>& prdTruth);

    virtual void buildDetailedTrackTruth(std::vector<DetailedTrackTruth> *output,
                        const Muon::MuonPatternCombination &pattern,
                        const std::vector<const PRD_MultiTruthCollection*>& prdTruth);

    void buildDetailedTrackTruthFromSegments(std::vector<DetailedSegmentTruth> *output,
                                             const Muon::MuonSegment &segment,
                                             const std::vector<const PRD_MultiTruthCollection*>& prdTruth);

  private:
    typedef InverseMultiMap<PRD_MultiTruthCollection> PRD_InverseTruth;

    ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc {this, "MuonIdHelperSvc", "Muon::MuonIdHelperSvc/MuonIdHelperSvc"};

    ToolHandle<Trk::ITruthTrajectoryBuilder> m_truthTrackBuilder{this,"TruthTrajectoryTool","Trk::ElasticTruthTrajectoryBuilder"};
    ToolHandle<Muon::IMdtDriftCircleOnTrackCreator> m_mdtCreator{this,"MdtDriftCircleCreator","Muon::MdtDriftCircleOnTrackCreator/MdtDriftCircleOnTrackCreator"};
    ToolHandle<Muon::IMuonClusterOnTrackCreator> m_muonClusterCreator{this,"MuonClusterCreator","Muon::MuonClusterOnTrackCreator/MuonClusterOnTrackCreator"};

    SubDetHitStatistics::SubDetType findSubDetType(Identifier id);

    const MuonSimData::Deposit* getDeposit(const MuonSimDataCollection& simCol, 
                                                                            HepMC::ConstGenParticlePtr genPart, 
                                                                            const Identifier& id);
    const MuonSimDataCollection* retrieveTruthCollection( std::string colName );
 
    void addTrack(DetailedMuonPatternTruthCollection *output,
		  const ElementLink<DataVector<Muon::MuonPatternCombination> > &track,
		  const std::vector<const PRD_MultiTruthCollection*>& orderedPRD_Truth,
		  const PRD_InverseTruth& inverseTruth);

    void addDetailedTrackTruth(std::vector<DetailedTrackTruth> *output,
                               const Muon::MuonPatternCombination& pattern,
                               const std::vector<const PRD_MultiTruthCollection*>& orderedPRD_Truth,
                               const PRD_InverseTruth& inverseTruth);

    SubDetHitStatistics countPRDsOnTruth(const TruthTrajectory& traj,
					 const PRD_InverseTruth& inverseTruth,
                                         std::set<Muon::MuonStationIndex::ChIndex> chIndices);

    Amg::Vector3D getPRDTruthPosition(const Muon::MuonSegment &segment, 
				      std::list<HepMC::ConstGenParticlePtr> genPartList,
				      int truthPos,
				      std::set<Muon::MuonStationIndex::ChIndex> chIndices);

    void addDetailedTrackTruthFromSegment(std::vector<DetailedSegmentTruth> *output,
                                          const Muon::MuonSegment &segment,
                                          const std::vector<const PRD_MultiTruthCollection*>& orderedPRD_Truth,
                                          const PRD_InverseTruth& inverseTruth);

  };
  
} // end namespace Trk

#endif
