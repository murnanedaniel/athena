/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MOOSEGMENTFINDERS_MUOSEGMENTFINDERALGS_H
#define MOOSEGMENTFINDERS_MUOSEGMENTFINDERALGS_H

#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "CscSegmentMakers/ICscSegmentFinder.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "MuonIdHelpers/IMuonIdHelperSvc.h"
#include "MuonPattern/MuonPatternChamberIntersect.h"
#include "MuonPrepRawData/CscPrepDataCollection.h"
#include "MuonPrepRawData/MdtPrepDataCollection.h"
#include "MuonPrepRawData/RpcPrepDataCollection.h"
#include "MuonPrepRawData/TgcPrepDataCollection.h"
#include "MuonRecHelperTools/MuonEDMPrinterTool.h"
#include "MuonRecToolInterfaces/IMuonClusterOnTrackCreator.h"
#include "MuonRecToolInterfaces/IMuonSegmentMaker.h"
#include "MuonSegment/MuonSegmentCombinationCollection.h"
#include "MuonSegmentMakerToolInterfaces/IMuonClusterSegmentFinder.h"
#include "MuonSegmentMakerToolInterfaces/IMuonClusterSegmentFinderTool.h"
#include "MuonSegmentMakerToolInterfaces/IMuonPatternCalibration.h"
#include "MuonSegmentMakerToolInterfaces/IMuonPatternSegmentMaker.h"
#include "MuonSegmentMakerToolInterfaces/IMuonSegmentOverlapRemovalTool.h"
#include "TrkSegment/SegmentCollection.h"
#include "TrkTruthData/PRD_MultiTruthCollection.h"

class MuonSegmentFinderAlg : public AthReentrantAlgorithm {
public:
    MuonSegmentFinderAlg(const std::string& name, ISvcLocator* pSvcLocator);

    virtual ~MuonSegmentFinderAlg() = default;

    virtual StatusCode initialize() override;
    virtual StatusCode execute(const EventContext& ctx) const override;

private:
    ServiceHandle<Muon::IMuonIdHelperSvc> m_idHelperSvc{
        this,
        "MuonIdHelperSvc",
        "Muon::MuonIdHelperSvc/MuonIdHelperSvc",
    };

    ToolHandle<Muon::MuonEDMPrinterTool> m_printer{
        this,
        "EDMPrinter",
        "Muon::MuonEDMPrinterTool/MuonEDMPrinterTool",
    };  //<! helper printer tool
    ToolHandle<Muon::IMuonPatternCalibration> m_patternCalibration{
        this,
        "MuonPatternCalibration",
        "Muon::MuonPatternCalibration/MuonPatternCalibration",
    };
    ToolHandle<Muon::IMuonPatternSegmentMaker> m_patternSegmentMaker{
        this,
        "MuonPatternSegmentMaker",
        "Muon::MuonPatternSegmentMaker/MuonPatternSegmentMaker",
    };
    ToolHandle<Muon::IMuonSegmentMaker> m_segmentMaker{
        this,
        "SegmentMaker",
        "Muon::DCMathSegmentMaker/DCMathSegmentMaker",
    };
    ToolHandle<Muon::IMuonClusterSegmentFinder> m_clusterSegMaker{
        this,
        "MuonClusterSegmentFinder",
        "Muon::MuonClusterSegmentFinder/MuonClusterSegmentFinder",
    };
    ToolHandle<Muon::IMuonSegmentOverlapRemovalTool> m_segmentOverlapRemovalTool{
        this,
        "MuonSegmentOverlapRemovalTool",
        "Muon::MuonSegmentOverlapRemovalTool/MuonSegmentOverlapRemovalTool",
    };
    ToolHandle<Muon::IMuonClusterOnTrackCreator> m_clusterCreator{
        this,
        "ClusterCreator",
        "Muon::MuonClusterOnTrackCreator/MuonClusterOnTrackCreator",
    };  //<! pointer to muon cluster rio ontrack creator
    ToolHandle<Muon::IMuonClusterOnTrackCreator> m_mmClusterCreator{
        this,
        "MMClusterCreator",
        "",
    };  //<! pointer to mm cluster rio ontrack creator
    ToolHandle<Muon::IMuonClusterSegmentFinderTool> m_clusterSegMakerNSW{
        this,
        "MuonClusterSegmentFinderTool",
        "",
    };
    ToolHandle<ICscSegmentFinder> m_csc2dSegmentFinder{
        this,
        "Csc2dSegmentMaker",
        "Csc2dSegmentMaker/Csc2dSegmentMaker",
    };
    ToolHandle<ICscSegmentFinder> m_csc4dSegmentFinder{
        this,
        "Csc4dSegmentMaker",
        "Csc4dSegmentMaker/Csc4dSegmentMaker",
    };

    // the following Trk::SegmentCollection MuonSegments are standard MuonSegments, the MuGirl segments are stored in MuonCreatorAlg.h
    SG::WriteHandleKey<Trk::SegmentCollection> m_segmentCollectionKey{
        this,
        "SegmentCollectionName",
        "TrackMuonSegments",
        "Muon Segments",
    };
    SG::WriteHandleKey<Trk::SegmentCollection> m_segmentNSWCollectionKey{ //this collection of segments are used to perform the alignment of the NSW
      this,
        "NSWSegmentCollectionName",
        "TrackMuonNSWSegments",
        "WriteHandleKey for NSW Segments",
    };
    SG::ReadHandleKey<Muon::CscPrepDataContainer> m_cscPrdsKey{
        this,
        "CSC_clusterkey",
        "CSC_Clusters",
        "CSC PRDs",
    };
    SG::ReadHandleKey<Muon::MdtPrepDataContainer> m_mdtPrdsKey{
        this,
        "MDT_PRDs",
        "MDT_DriftCircles",
        "MDT PRDs",
    };
    SG::ReadHandleKey<Muon::RpcPrepDataContainer> m_rpcPrdsKey{
        this,
        "RPC_PRDs",
        "RPC_Measurements",
        "RPC PRDs",
    };
    SG::ReadHandleKey<Muon::TgcPrepDataContainer> m_tgcPrdsKey{
        this,
        "TGC_PRDs",
        "TGC_Measurements",
        "TGC PRDs",
    };
    SG::ReadHandleKey<MuonPatternCombinationCollection> m_patternCollKey{
        this,
        "MuonLayerHoughCombisKey",
        "MuonLayerHoughCombis",
        "Hough combinations",
    };
    SG::ReadHandleKey<PRD_MultiTruthCollection> m_tgcTruth{
        this,
        "TGCTruth",
        "TGC_TruthMap",
        "TGC PRD Multi-truth Collection",
    };
    SG::ReadHandleKey<PRD_MultiTruthCollection> m_rpcTruth{
        this,
        "RPCTruth",
        "RPC_TruthMap",
        "RPC PRD Multi-truth Collection",
    };

    void createSegmentsWithMDTs(const Muon::MuonPatternCombination* patt, Trk::SegmentCollection* segs,
                                const std::vector<const Muon::RpcPrepDataCollection*>& rpcCols,
                                const std::vector<const Muon::TgcPrepDataCollection*>& tgcCols, const EventContext& ctx) const;
    void createSegmentsFromClusters(const EventContext& ctx, const Muon::MuonPatternCombination* patt, Trk::SegmentCollection* segments, Trk::SegmentCollection* segmentsNSW) const;

    Gaudi::Property<bool> m_printSummary{this, "PrintSummary", false};
    Gaudi::Property<bool> m_doTGCClust{this, "doTGCClust", false, "selection flags for cluster based segment finding"};
    Gaudi::Property<bool> m_doRPCClust{this, "doRPCClust", false, "selection flags for cluster based segment finding"};
    Gaudi::Property<bool> m_doClusterTruth{this, "doClusterTruth", false, "selection flags for cluster based segment finding"};
};

#endif
