/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONTRACKSTATISTICSALG_MUONTRACKSTATISTICSALG_H
#define MUONTRACKSTATISTICSALG_MUONTRACKSTATISTICSALG_H

#include <fstream>
#include <string>

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "StoreGate/ReadHandleKeyArray.h"
#include "TrkTrack/TrackCollection.h"
#include "TrkTruthData/DetailedTrackTruthCollection.h"

class MuonTrackStatisticsTool;

class MuonTrackStatisticsAlg : public AthAlgorithm {
public:
    // Algorithm Constructor
    MuonTrackStatisticsAlg(const std::string &name, ISvcLocator *pSvcLocator);
    ~MuonTrackStatisticsAlg(){};

    // Gaudi algorithm hooks
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:
    // ServiceHandle<Muon::IMuonEDMHelperSvc>   m_edmHelperSvc;
    ToolHandle<MuonTrackStatisticsTool> m_statisticsTool;

    /** name of external file to write statistics */
    std::string m_fileName;

    /** output file*/
    std::ofstream m_fileOutput;

    bool m_writeToFile;
    bool m_doTruth;

    SG::ReadHandleKeyArray<TrackCollection> m_trackKeys{
        this, "TrackLocationList", {"MuonSpectrometerTracks"}, "track collections to track"};
    SG::ReadHandleKeyArray<DetailedTrackTruthCollection> m_truthKeys{
        this, "TruthTrackLocationList", {"MuonSpectrometerTracksTruth"}, "truth track collections"};

    void storeTruthTracks(void);
};

#endif  // MUONTRACKSTATISTICSALG_MUONTRACKSTATISTICSALG_H
