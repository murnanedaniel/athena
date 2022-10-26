/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "MdtMathSegmentFinder.h"

#include "MuonReadoutGeometry/MdtReadoutElement.h"
#include "TrkDriftCircleMath/DCSLFitter.h"
#include "TrkDriftCircleMath/DCStatistics.h"
#include "TrkDriftCircleMath/Road.h"
#include "TrkDriftCircleMath/SegmentFinder.h"

 

namespace Muon {

    MdtMathSegmentFinder::MdtMathSegmentFinder(const std::string& t, const std::string& n, const IInterface* p) : AthAlgTool(t, n, p) {
        declareInterface<IMdtSegmentFinder>(this);

        declareProperty("FinderDebugLevel", m_finderDebugLevel = 0, "switch on debug output of finder");
        declareProperty("DoDrop", m_doDrop = true, "Recursive outlier removal");
        declareProperty("Chi2PerDofDropping", m_chi2PerDofDrop = 10., "Chi2 cut for recursive outlier removal");
        declareProperty("RatioEmptyTubeCut", m_ratioEmptyTubesCut = 1.1, "holes/hits cut - holes are all non-hits along the line");
        declareProperty("EnableSeedCleaning", m_enableSeedCleaning = false, "only use dc witout neighbours as seeds");
        declareProperty("OccupancyThreshold", m_occupancyThreshold = 0.3, "occupancy threshold before enabling seed cleaning");
        declareProperty("OccupancyCutoff", m_occupancyCutOff = 0.8, "above the occupancy threshold no segment finding");
        declareProperty("RPCAssocationPullcut", m_rpcAssociationPullCut = 5., "Association cut for RPCs");
        declareProperty("TGCAssocationPullcut", m_tgcAssociationPullCut = 5., "Association cut for TGCs");
        declareProperty("MDTAssocationPullcut", m_mdtAssociationPullCut = 5., "Association cut for MDTs");
        declareProperty("UseChamberTheta", m_useChamberTheta = true, "Always look for pointing segments");
        declareProperty("SortSegmentWithAllHits", m_doAllHitSort = true, "Including triggers in segment selection");
        declareProperty("AssociationRoadWidth", m_roadWidth = 1.5, "Road width used during hit association with seed lines");
        declareProperty("SearchForPointingSegments", m_useChamberTheta = true, "Always search for pointing segments");
        declareProperty("DoRoadSeeding", m_doRoadAngleSeeding = true);
        declareProperty("DoIPSeeding", m_doIPAngleSeeding = true);
        declareProperty("TightRoadCut", m_tightRoadCut = 0.1);
        declareProperty("MaxHitsPerFullSearch", m_maxHitsPerFullSearch = 100);
        declareProperty("DoSingleMultiLayerScan", m_doSingleMultiLayerScan = true, "Look for segments in one multi layer");
        declareProperty("RecoverMdtOutliers", m_recoverMdtOutliers = true, "Recover MDT outliers after fit");
        declareProperty("RemoveSingleMdtOutliers", m_removeSingleOutliers = true, "Remove single MDT outliers");
        declareProperty("DoCurvedSegmentFinder", m_doCurvedSegmentFinder = false, "Use the curved segment finding routine");
        declareProperty("DeltaCutT0Segments", m_deltaCutT0Segments = 5., "Delta cut for segments with T0 fit");
        declareProperty("ResidualCutT0Segments", m_residualCutT0Segments = 1., "Residual cut for segments with T0 fit");
        declareProperty("UseSegmentQuality", m_useSegmentQuality = false, "Use segment quality in hit dropping");
    }

    StatusCode MdtMathSegmentFinder::initialize() {
        if (!m_dcslFitProvider.empty()) {
            ATH_CHECK(m_dcslFitProvider.retrieve());
            ATH_MSG_INFO(" Using fitter from " << m_dcslFitProvider);
        }
        ATH_CHECK(m_idHelperSvc.retrieve());
        return StatusCode::SUCCESS;
    }

    const TrkDriftCircleMath::SegVec MdtMathSegmentFinder::findSegments(const TrkDriftCircleMath::DCVec& dcvec,
                                                                        const TrkDriftCircleMath::CLVec& clvec,
                                                                        const TrkDriftCircleMath::Road& road,
                                                                        const TrkDriftCircleMath::DCStatistics& dcstats,
                                                                        const TrkDriftCircleMath::ChamberGeometry* multiGeo = nullptr) const {
        // setup finder
        std::unique_ptr<TrkDriftCircleMath::SegmentFinder> segmentFinder = std::make_unique<TrkDriftCircleMath::SegmentFinder>(m_roadWidth, m_mdtAssociationPullCut, false);

        // set debug level
        segmentFinder->debugLevel(m_finderDebugLevel);

        // configure uasge of chamber position for angular seeding
        segmentFinder->setUseChamberPhi(m_useChamberTheta);

        // enable dropping of hits
        segmentFinder->setDropHits(m_doDrop);

        // enable seed cleaing
        segmentFinder->setSeedCleaning(m_enableSeedCleaning);

        // do single multilayer scan?
        segmentFinder->setSingleMultiLayerScan(m_doSingleMultiLayerScan);

        // set chi2/ndof threshold for cleaning of segments
        segmentFinder->setChi2DropCut(m_chi2PerDofDrop);

        // set ratio for dropping segments with many empty tubes
        segmentFinder->setRatioEmptyTubesCut(m_ratioEmptyTubesCut);

        // set sort mode segment finder
        segmentFinder->setSortSegmentsUsingAllHits(m_doAllHitSort);

        // set RPC pull cut
        segmentFinder->setRPCPullCut(m_rpcAssociationPullCut);

        // set TGC pull cut
        segmentFinder->setTGCPullCut(m_tgcAssociationPullCut);

        // set MDT outlier recovery
        segmentFinder->setRecoverMDT(m_recoverMdtOutliers);

        // set removal of single outliers
        segmentFinder->setRemoveSingleOutliers(m_removeSingleOutliers);

        // set the curved segment finder
        segmentFinder->setCurvedSegmentFinder(m_doCurvedSegmentFinder);

        // set removal of single outliers
        segmentFinder->setDeltaCutT0(m_deltaCutT0Segments);  // reset defaults

        // set removal of single outliers
        segmentFinder->setResidualCutT0(m_residualCutT0Segments);

        // set use of segment quality
        segmentFinder->setUseSegmentQuality(m_useSegmentQuality);

        if (!m_dcslFitProvider.empty()) {
            std::shared_ptr<const TrkDriftCircleMath::DCSLFitter> fitter(m_dcslFitProvider->getFitter(), Muon::IDCSLFitProvider::Unowned{});
            segmentFinder->setFitter(fitter);
        } else {
            segmentFinder->setFitter(std::make_shared<TrkDriftCircleMath::DCSLFitter>());
        }

        // set angle prediction from road
        segmentFinder->setPhiRoad(road.angle(), road.chamberAngle(), road.width(), m_doRoadAngleSeeding, m_doIPAngleSeeding);

        // set pointer to geometry
        segmentFinder->setMdtGeometry(multiGeo);

        // set seed cleaning
        bool highOccupancy = false;
        bool aboveOccupancyCut = false;
        double occupancyMax = 0.;
        unsigned int nmdtHits = 0;
        // calculate multi layer occupancy
        TrkDriftCircleMath::DCStatCit mit = dcstats.begin();
        TrkDriftCircleMath::DCStatCit mit_end = dcstats.end();
        for (; mit != mit_end; ++mit) {
            unsigned int channels = mit->first->getNLayers() * mit->first->getNtubesperlayer();
            double occupancy = (double)mit->second / (double)channels;

            nmdtHits += mit->second;

            if (occupancy > occupancyMax) occupancyMax = occupancy;

            if (occupancy > m_occupancyThreshold) highOccupancy = true;

            if (occupancy > m_occupancyCutOff) aboveOccupancyCut = true;
        }

        // sanity check
        if (nmdtHits != dcvec.size()) {
            ATH_MSG_WARNING(" inconsistent number of mdt hits " << nmdtHits << " from vec " << dcvec.size());
            nmdtHits = dcvec.size();
        }

        if (aboveOccupancyCut) {
            ATH_MSG_DEBUG(" layer with occupancy above threshold, aborting segment finding "
                          << occupancyMax << " cut " << m_occupancyCutOff << " nhits " << nmdtHits << " cut " << m_maxHitsPerFullSearch);
            return TrkDriftCircleMath::SegVec();
        }

        // enable seed cleaning
        if (highOccupancy || nmdtHits > m_maxHitsPerFullSearch) {
            ATH_MSG_DEBUG(" switch to fast search: occupancy " << occupancyMax << " cut " << m_occupancyThreshold << " nhits " << nmdtHits
                                                               << " cut " << m_maxHitsPerFullSearch);

            // to speed up reconstruction use default fitter
            if (!m_dcslFitProvider.empty()) {
                segmentFinder->setFitter(std::make_shared<TrkDriftCircleMath::DCSLFitter>());
            }

            // use tight road cuts and only look for pointing segments
            segmentFinder->setPhiRoad(road.chamberAngle(), road.chamberAngle(), m_tightRoadCut);
        }

        return segmentFinder->findSegments(dcvec, clvec);
    }

}  // namespace Muon
