/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "DCMathSegmentMaker.h"

#include <cassert>
#include <functional>
#include <iostream>

#include "AthenaKernel/Timeout.h"
#include "EventPrimitives/EventPrimitivesHelpers.h"
#include "EventPrimitives/EventPrimitivesToStringConverter.h"
#include "FourMomUtils/xAODP4Helpers.h"
#include "GeoPrimitives/GeoPrimitivesToStringConverter.h"
#include "MuonCompetingRIOsOnTrack/CompetingMuonClustersOnTrack.h"
#include "MuonPrepRawData/MdtPrepData.h"
#include "MuonPrepRawData/RpcPrepData.h"
#include "MuonPrepRawData/TgcPrepData.h"
#include "MuonRIO_OnTrack/MdtDriftCircleOnTrack.h"
#include "MuonRIO_OnTrack/RpcClusterOnTrack.h"
#include "MuonRIO_OnTrack/TgcClusterOnTrack.h"
#include "MuonReadoutGeometry/MdtReadoutElement.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/RpcReadoutElement.h"
#include "MuonSegment/MuonSegmentQuality.h"
#include "TrkDriftCircleMath/ClusterId.h"
#include "TrkDriftCircleMath/DriftCircle.h"
#include "TrkDriftCircleMath/MdtMultiChamberGeometry.h"
#include "TrkDriftCircleMath/ResidualWithSegment.h"
#include "TrkDriftCircleMath/Road.h"
#include "TrkDriftCircleMath/Segment.h"
#include "TrkDriftCircleMath/SegmentFinder.h"
#include "TrkEventPrimitives/DefinedParameter.h"
#include "TrkEventPrimitives/FitQuality.h"
#include "TrkEventPrimitives/JacobianCotThetaPtToThetaP.h"
#include "TrkEventPrimitives/LocalDirection.h"
#include "TrkEventPrimitives/LocalParameters.h"
#include "TrkEventPrimitives/ParamDefs.h"
#include "TrkParameters/TrackParameters.h"
#include "TrkRIO_OnTrack/RIO_OnTrack.h"
#include "TrkRoad/TrackRoad.h"
#include "TrkSurfaces/PlaneSurface.h"
#include "TrkSurfaces/StraightLineSurface.h"
#include "TrkTrack/Track.h"
namespace {

    double cot(double x) {
        /// use a cut off scale between
        constexpr double dX = std::numeric_limits<float>::epsilon();
        if (std::abs(x) < dX || std::abs(x - M_PI) < dX) return std::numeric_limits<float>::max();
        /// Use simple relation between cot(x) = tan( PI/2 -x) to convert into each other
        return std::tan(M_PI_2 - x);
    }

}  // namespace
namespace Muon {

    DCMathSegmentMaker::DCMathSegmentMaker(const std::string& t, const std::string& n, const IInterface* p) :
        AthAlgTool(t, n, p)  {
        declareInterface<IMuonSegmentMaker>(this);
        declareInterface<IMuonSegmentTriggerHitAssociator>(this);       
    }

    StatusCode DCMathSegmentMaker::initialize() {
        // retrieve MuonDetectorManager
        ATH_CHECK(m_DetectorManagerKey.initialize());
        ATH_CHECK(m_mdtCreator.retrieve());
        ATH_CHECK(m_mdtCreatorT0.retrieve());
        ATH_CHECK(m_clusterCreator.retrieve());
        ATH_CHECK(m_compClusterCreator.retrieve());
        ATH_CHECK(m_idHelperSvc.retrieve());
        ATH_CHECK(m_printer.retrieve());
        ATH_CHECK(m_edmHelperSvc.retrieve());
        ATH_CHECK(m_segmentFinder.retrieve());
        ATH_CHECK(m_segmentSelectionTool.retrieve());

        if (m_refitParameters) { ATH_CHECK(m_segmentFitter.retrieve()); }

        if (!m_dcslFitProvider.empty()) {
            ATH_CHECK(m_dcslFitProvider.retrieve());
            ATH_MSG_INFO(" Using " << m_dcslFitProvider);
        }

        // initialise for data handles
        ATH_CHECK(m_rpcKey.initialize());
        ATH_CHECK(m_tgcKey.initialize());
        ATH_CHECK(m_mdtKey.initialize());
        ATH_CHECK(m_chamberGeoKey.initialize());
       
        return StatusCode::SUCCESS;
    }

    MuonSegment* DCMathSegmentMaker::associateTriggerHits(const MuonSegment& seg, const std::vector<const MuonClusterOnTrack*>& clus,
                                                                bool includeEtaHits) const {
        const EventContext& ctx = Gaudi::Hive::currentContext();
        ATH_MSG_DEBUG("associateTriggerHits: clusters " << clus.size());
        // extract hits
        std::vector<const MdtDriftCircleOnTrack*> mdts;
        mdts.reserve(seg.containedMeasurements().size());
        for (const Trk::MeasurementBase* hit : seg.containedMeasurements()) {
            const MdtDriftCircleOnTrack* mdt = dynamic_cast<const MdtDriftCircleOnTrack*>(hit);
            if (mdt) mdts.push_back(mdt);
        }

        // global to local transformation for chamber
        Amg::Transform3D gToStation = seg.associatedSurface().transform().inverse();

        // set to get Identifiers of chambers with hits
        std::unique_ptr<TrkDriftCircleMath::MdtMultiChamberGeometry> multiGeo;
        if (m_doGeometry) {
            std::set<Identifier> chamberSet;
            for (const MdtDriftCircleOnTrack* drift_circle : mdts) {
                Identifier dcid = drift_circle->identify();
                Identifier elId = m_idHelperSvc->mdtIdHelper().elementID(dcid);
                chamberSet.insert(elId);
            }

            // vector to store chamber geometries
            std::vector<TrkDriftCircleMath::MdtChamberGeometry> geos;

            // loop over chambers
            geos.reserve(chamberSet.size());
        for (const Identifier& cham_id : chamberSet) { geos.push_back(createChamberGeometry(cham_id, gToStation)); }

            // create new geometry
            multiGeo = std::make_unique<TrkDriftCircleMath::MdtMultiChamberGeometry>(geos);
        }

        if (mdts.size() < 2) return nullptr;
        const MdtDriftCircleOnTrack* firstRot = findFirstRotInChamberWithMostHits(mdts);

        if (!firstRot) { return nullptr; }

        const MuonGM::MdtReadoutElement* detEl = firstRot->detectorElement();

        if (!detEl) {
            ATH_MSG_WARNING(" no MdtReadoutElement found, returning 0 ");
            return nullptr;
        }

        // identifier
        Identifier chid = firstRot->identify();

        // endcap or barrel
        bool isEndcap = m_idHelperSvc->mdtIdHelper().isEndcap(chid);

        // define axis of chamber in global coordinates
        Amg::Transform3D amdbToGlobal = seg.associatedSurface().transform();

        std::vector<const MuonClusterOnTrack*> clusters = clus;
        std::vector<std::unique_ptr<const Trk::MeasurementBase> > measToBeDeleted;
        // if requested, retrieve trigger eta hits
        if (includeEtaHits) {
            addEtaHits(clusters, measToBeDeleted, isEndcap);
            ATH_MSG_DEBUG("adding eta hits: original clusters " << clus.size() << " extended " << clusters.size());
        }

        /* ***** create cluster hits ******** */
        ClusterVecPair spVecs;
        if (m_doSpacePoints)
            spVecs = create2DClusters(clusters);
        else
            spVecs = create1DClusters(clusters);
        TrkDriftCircleMath::CLVec cls = createClusterVec(chid, spVecs.first, gToStation);

        /* ***** create MDT hits ************ */
        ATH_MSG_DEBUG(" adding mdts " << mdts.size());

        TrkDriftCircleMath::DCOnTrackVec dcs;
        // reset chamber statistics
        bool firstMdt = true;
        double phimin{9999.}, phimax{-9999.};
        unsigned index = 0;
        for (const Trk::MeasurementBase* hit : seg.containedMeasurements()) {
            const MdtDriftCircleOnTrack* mdt = dynamic_cast<const MdtDriftCircleOnTrack*>(hit);

            if (!mdt) { continue; }
            Identifier id = mdt->identify();

            // calculate local AMDB position
            Amg::Vector3D locPos = gToStation * mdt->prepRawData()->globalPosition();
            TrkDriftCircleMath::LocVec2D lpos(locPos.y(), locPos.z());

            double r = mdt->localParameters()[Trk::locR];
            double dr = Amg::error(mdt->localCovariance(), Trk::locR);

            // create identifier
            TrkDriftCircleMath::MdtId mdtid(m_idHelperSvc->mdtIdHelper().isBarrel(id), m_idHelperSvc->mdtIdHelper().multilayer(id) - 1,
                                            m_idHelperSvc->mdtIdHelper().tubeLayer(id) - 1, m_idHelperSvc->mdtIdHelper().tube(id) - 1);

            // create new DriftCircle
            TrkDriftCircleMath::DriftCircle dc(lpos, r, dr, TrkDriftCircleMath::DriftCircle::InTime, mdtid, index, mdt);
            TrkDriftCircleMath::DCOnTrack dcOnTrack(dc, 1., 1.);
            ATH_MSG_VERBOSE(" new MDT hit " << m_idHelperSvc->toString(id) << " x " << lpos.x() << " y " << lpos.y());

            TubeEnds tubeEnds = localTubeEnds(*mdt, gToStation, amdbToGlobal);
            if (firstMdt) {
                phimin = tubeEnds.phimin;
                phimax = tubeEnds.phimax;
                firstMdt = false;
            } else {
                updatePhiRanges(tubeEnds.phimin, tubeEnds.phimax, phimin, phimax);
            }

            dcs.push_back(dcOnTrack);
            ++index;
        }

        // create line
        double angleYZ = seg.localDirection().angleYZ();
        const Amg::Vector3D lpos = gToStation * seg.globalPosition();
        TrkDriftCircleMath::LocVec2D segPos(lpos.y(), lpos.z());
        TrkDriftCircleMath::Line segPars(segPos, angleYZ);

        const Trk::FitQuality* fq = seg.fitQuality();
        double chi2 = fq ? fq->chiSquared() : 0.;
        int ndof = fq ? fq->numberDoF() : 1;

        TrkDriftCircleMath::Segment segment(segPars, dcs, chi2, ndof, Amg::error(seg.localCovariance(), Trk::locY),
                                            Amg::error(seg.localCovariance(), Trk::theta));
        segment.clusters(cls);
        segment.hitsOnTrack(6);  // hack just putting a number above the cut of at least 2

        segmentCreationInfo sInfo(spVecs, multiGeo.get(), gToStation, amdbToGlobal, phimin, phimax);
        std::unique_ptr<MuonSegment> newSeg = createSegment(ctx, segment, chid, seg.globalPosition(), seg.globalDirection(), mdts, true, sInfo);

        return newSeg.release();
    }

    void DCMathSegmentMaker::addEtaHits(std::vector<const MuonClusterOnTrack*>& clusters,
                                        std::vector<std::unique_ptr<const Trk::MeasurementBase> >& measurementsToBeDeleted,
                                        bool isEndcap) const {
        std::set<IdentifierHash> chIdHs;
        for (const MuonClusterOnTrack* clus : clusters) chIdHs.insert(clus->collectionHash());
        ATH_MSG_VERBOSE("Number of collections " << chIdHs.size());

        if (isEndcap) {
            const Muon::TgcPrepDataContainer* prdContainer = nullptr;
            SG::ReadHandle<Muon::TgcPrepDataContainer> TgcCont(m_tgcKey);
            if (!TgcCont.isValid()) {
                ATH_MSG_WARNING("Cannot retrieve TgcPrepDataContainer ");
                return;
            }
            prdContainer = TgcCont.cptr();
            if (!prdContainer) {
                ATH_MSG_WARNING("No TGC prd container retrieved");
                return;
            }

            // loop over chambers and get collections
            for (const IdentifierHash& chit : chIdHs) {
                auto collptr = prdContainer->indexFindPtr(chit);
                if (!collptr) {
                    ATH_MSG_VERBOSE("Could not find collection ");
                    continue;
                }
                const TgcPrepDataCollection& col = *collptr;
                ATH_MSG_VERBOSE("Found collection: " << m_idHelperSvc->toStringChamber(col.identify()) << " size " << col.size());
                TgcPrepDataCollection::const_iterator hit = col.begin();
                TgcPrepDataCollection::const_iterator hit_end = col.end();
                for (; hit != hit_end; ++hit) {
                    // skip phi hits
                    if (m_idHelperSvc->measuresPhi((*hit)->identify())) continue;

                    const MuonClusterOnTrack* clus = m_clusterCreator->createRIO_OnTrack(**hit, (*hit)->globalPosition());
                    if (clus) {
                        ATH_MSG_VERBOSE("  adding hit: " << m_idHelperSvc->toString(clus->identify()));
                        clusters.push_back(clus);
                        measurementsToBeDeleted.emplace_back(clus);
                    }
                }
            }

        } else {
            const Muon::RpcPrepDataContainer* prdContainer = nullptr;
            SG::ReadHandle<Muon::RpcPrepDataContainer> RpcCont(m_rpcKey);
            if (!RpcCont.isValid()) {
                ATH_MSG_WARNING("Cannot retrieve RpcPrepDataContainer ");
                return;
            }
            prdContainer = RpcCont.cptr();
            if (!prdContainer) {
                ATH_MSG_WARNING("No RPC prd container retrieved");
                return;
            }

            // loop over chambers and get collections
            for (const IdentifierHash& chit : chIdHs) {
                auto collptr = prdContainer->indexFindPtr(chit);
                if (!collptr) continue;

                const RpcPrepDataCollection& col = *collptr;
                RpcPrepDataCollection::const_iterator hit = col.begin();
                RpcPrepDataCollection::const_iterator hit_end = col.end();
                for (; hit != hit_end; ++hit) {
                    // skip phi hits
                    if (m_idHelperSvc->measuresPhi((*hit)->identify())) continue;
                    const MuonClusterOnTrack* clus = m_clusterCreator->createRIO_OnTrack(**hit, (*hit)->globalPosition());
                    if (clus) {
                        clusters.push_back(clus);
                        measurementsToBeDeleted.emplace_back(clus);
                    }
                }
            }
        }
    }

    void DCMathSegmentMaker::find(const Amg::Vector3D& roadpos, const Amg::Vector3D& roaddir,
                                  const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                                  const std::vector<const MuonClusterOnTrack*>& clusters, bool hasPhiMeasurements,
                                  Trk::SegmentCollection* segColl, double momentum, double sinAngleCut) const {
        const EventContext& ctx = Gaudi::Hive::currentContext();
        if (m_doTimeOutChecks && Athena::Timeout::instance().reached()) {
            ATH_MSG_DEBUG("Timeout reached. Aborting sequence.");
            return;
        }

        ATH_MSG_DEBUG("In find, passed " << mdts.size() << " RIO_OnTracks");

        if (mdts.size() < 3) return;

        const MdtDriftCircleOnTrack* firstRot = findFirstRotInChamberWithMostHits(mdts);

        if (!firstRot) { return; }

        const MuonGM::MdtReadoutElement* detEl = firstRot->detectorElement();

        if (!detEl) {
            ATH_MSG_WARNING(" no MdtReadoutElement found, returning 0 ");
            return;
        }

        // identifier
        Identifier chid = firstRot->identify();

        // endcap or barrel
        bool isEndcap = m_idHelperSvc->mdtIdHelper().isEndcap(chid);

        // global to local transformation for chamber
        Amg::Transform3D gToStation = detEl->GlobalToAmdbLRSTransform();

        // define axis of chamber in global coordinates
        Amg::Transform3D amdbToGlobal = detEl->AmdbLRSToGlobalTransform();

        // transform nominal pointing chamber position into surface frame
        Amg::Vector3D globalDirCh = Amg::Vector3D(detEl->center().x(), detEl->center().y(), detEl->center().z());
        Amg::Vector3D dirCh(gToStation.linear() * globalDirCh);
        double chamber_angleYZ = atan2(dirCh.z(), dirCh.y());

        Amg::Vector3D roaddir2 = roaddir;
        double dotprod = globalDirCh.perp() * std::sin(roaddir2.theta()) + globalDirCh.z() * std::cos(roaddir2.theta());
        if (dotprod < 0) roaddir2 = -roaddir2;

        // transform the global direction into the surface frame
        Amg::Vector3D d(gToStation.linear() * roaddir2);
        // calculate the local road angles in the surface frame
        double road_angleXZ = std::atan2(d.z(), d.x());
        double road_angleYZ = std::atan2(d.z(), d.y());

        if (!hasPhiMeasurements) road_angleXZ = M_PI;  // if no phi, take phi perpendicular to plane
        ATH_MSG_VERBOSE("global road dir " << Amg::toString(roaddir2) << " XZ " << road_angleXZ << " YZ " << road_angleYZ << " isEndcap "
                                           << isEndcap << " central phi " << detEl->center().phi() << " r " << detEl->center().perp()
                                           << " z " << detEl->center().z());

        // rescale errors for low momentum
        double errorScale = errorScaleFactor(chid, momentum, hasPhiMeasurements);

        /* ***** create cluster hits ******** */
        ATH_MSG_DEBUG(" adding clusters " << clusters.size());
        ClusterVecPair spVecs;
        if (m_doSpacePoints)
            spVecs = create2DClusters(clusters);
        else
            spVecs = create1DClusters(clusters);
        TrkDriftCircleMath::CLVec cls = createClusterVec(chid, spVecs.first, gToStation);

        /* ***** create MDT hits ************ */
        ATH_MSG_DEBUG(" adding mdts " << mdts.size());
        for (auto it : mdts) ATH_MSG_DEBUG(*it);

        // set to get Identifiers of chambers with hits
        std::set<Identifier> chamberSet;
        double phimin{-9999}, phimax{9999};
        TrkDriftCircleMath::DCStatistics dcStatistics;  // statistics on chamber occupancy
        TrkDriftCircleMath::DCVec dcs = createDCVec(mdts, errorScale, chamberSet, phimin, phimax, dcStatistics, gToStation, amdbToGlobal);

        // create geometry
        std::unique_ptr<TrkDriftCircleMath::MdtMultiChamberGeometry> multiGeo;
        if (m_doGeometry) {
            ATH_MSG_VERBOSE(" using chamber geometry with #chambers " << chamberSet.size());

            // vector to store chamber geometries
            std::vector<TrkDriftCircleMath::MdtChamberGeometry> geos;

            // loop over chambers
            geos.reserve(chamberSet.size());
            for (const Identifier& id : chamberSet) { geos.push_back(createChamberGeometry(id, gToStation)); }

            // create new geometry
            multiGeo = std::make_unique<TrkDriftCircleMath::MdtMultiChamberGeometry>(geos);
        }

        double angle = m_sinAngleCut;
        if (sinAngleCut > 0) angle = sinAngleCut;
        TrkDriftCircleMath::Road road(TrkDriftCircleMath::LocVec2D(0., 0.), road_angleYZ, chamber_angleYZ, angle);

        // call segment finder
        TrkDriftCircleMath::SegVec segs = m_segmentFinder->findSegments(dcs, cls, std::move(road), dcStatistics, multiGeo.get());

        ATH_MSG_DEBUG("Found " << segs.size() << " segments");

        // return
        if (segs.empty()) { return; }

        // loop over segments
        TrkDriftCircleMath::SegIt sit = segs.begin();
        TrkDriftCircleMath::SegIt sit_end = segs.end();
        segmentCreationInfo sInfo(spVecs, multiGeo.get(), gToStation, amdbToGlobal, phimin, phimax);
        for (; sit != sit_end; ++sit) {
            std::unique_ptr<MuonSegment> segment = createSegment(ctx, *sit, chid, roadpos, roaddir2, mdts, hasPhiMeasurements, sInfo);
            if (segment) segColl->push_back(segment.release());
        }
        ATH_MSG_DEBUG(" Done ");
    }

    std::unique_ptr<MuonSegment> DCMathSegmentMaker::createSegment(const EventContext& ctx, TrkDriftCircleMath::Segment& segment, const Identifier& chid,
                                                   const Amg::Vector3D& roadpos, const Amg::Vector3D& roaddir2,
                                                   const std::vector<const MdtDriftCircleOnTrack*>& mdts, bool hasPhiMeasurements,
                                                   segmentCreationInfo& sInfo) const {
        bool isEndcap = m_idHelperSvc->isEndcap(chid);
        // find all curved segments
        bool isCurvedSegment(false);
        MuonStationIndex::ChIndex chIndex = m_idHelperSvc->chamberIndex(chid);

        // MuonStationIndex::ChIndex chIndex;
        if (segment.hasCurvatureParameters()) {
            if (chIndex == Muon::MuonStationIndex::BIL || chIndex == MuonStationIndex::BML || chIndex == MuonStationIndex::BMS ||
                chIndex == MuonStationIndex::BOL)
                isCurvedSegment = true;
        }

        // remove segments with too few hits
        if (segment.hitsOnTrack() < 3) return nullptr;

        // convert segment parameters + x position from road
        const TrkDriftCircleMath::Line& line = segment.line();

        ATH_MSG_DEBUG("New segment: chi2 " << segment.chi2() << " ndof " << segment.ndof() << " line " << line.position().x() << ","
                                           << line.position().y() << " phi " << line.phi() << " associated clusters "
                                           << segment.clusters().size());

        // local position along x from road
        Amg::Vector3D lroadpos = sInfo.globalTrans * roadpos;
        Amg::Vector3D lroaddir = sInfo.globalTrans.linear() * roaddir2;

        // local x position of first tube used if no phi measurement is present
        double lxroad = 0.;

        if (hasPhiMeasurements) {
            // calculate local position of segment along tube using the road
            // calculate intersect pattern measurement plane
            double sphi = 0.;
            double cphi = lroaddir.x();
            // swap local y and z in the endcaps
            if (isEndcap) {
                sphi = lroaddir.y();
                lxroad = lroadpos.x() + (-lroadpos.y() + line.position().x()) * cphi / sphi;
            } else {
                sphi = lroaddir.z();
                lxroad = lroadpos.x() + (-lroadpos.z() + line.position().y()) * cphi / sphi;
            }

            double shortestTubeLen = 1e9;
            // loop over hits and get the shortest tube on the segment
            for (TrkDriftCircleMath::DCOnTrackCit dcit = segment.dcs().begin(); dcit != segment.dcs().end(); ++dcit) {
                if (dcit->state() != TrkDriftCircleMath::DCOnTrack::OnTrack) continue;

                const MdtDriftCircleOnTrack* riodc = mdts[dcit->index()];
                if (!riodc) continue;
                int lay = m_idHelperSvc->mdtIdHelper().tubeLayer(riodc->identify());
                int tube = m_idHelperSvc->mdtIdHelper().tube(riodc->identify());
                double tubelen = 0.5 * riodc->prepRawData()->detectorElement()->getActiveTubeLength(lay, tube);
                if (tubelen < shortestTubeLen) shortestTubeLen = tubelen;
            }
            // if the predicted position lies outside the chamber move it back inside
            if (std::abs(lxroad) > shortestTubeLen) {
                ATH_MSG_DEBUG("coordinates far outside chamber! using global position of first hit ");
                if (lxroad < 0.) shortestTubeLen *= -1.;
                lxroad = shortestTubeLen;
            }
        } else {
            lxroad = (sInfo.globalTrans * mdts[0]->prepRawData()->detectorElement()->surface(mdts[0]->identify()).center()).x();
        }

        // calculate local direction vector
        Amg::Vector3D lpos(lxroad, line.position().x(), line.position().y());

        // global position segment
        Amg::Vector3D gpos = sInfo.amdbTrans * lpos;

        // create new surface
        Amg::Transform3D surfaceTransform(sInfo.amdbTrans.rotation());
        surfaceTransform.pretranslate(gpos);
        double surfDim = 500.;
        std::unique_ptr<Trk::PlaneSurface> surf = std::make_unique<Trk::PlaneSurface>(surfaceTransform, surfDim, surfDim);

        // measurements
        Amg::Vector2D segLocPos(0., 0.);
        double linephi = line.phi();

        // now update the global direction using the local precision angle of the segment and the global phi angle of the
        // road.
        Amg::Vector3D gdir = updateDirection(linephi, *surf, roaddir2, isCurvedSegment);

        // extract RIO_OnTracks
        std::vector<std::pair<double,  std::unique_ptr<const Trk::MeasurementBase>> > rioDistVec;  // vector to store the distance of a ROT to the segment

        // associate MDT hits to segment
        std::set<Identifier> deltaVec;
        std::set<Identifier> outoftimeVec;
        
        associateMDTsToSegment(gdir, segment, mdts, sInfo.geom, sInfo.globalTrans, sInfo.amdbTrans, deltaVec, outoftimeVec, rioDistVec);
        std::vector<std::pair<double, std::unique_ptr<const Trk::MeasurementBase>>> garbage_collector;

        TrkDriftCircleMath::DCSLHitSelector hitSelector;

        if (m_redo2DFit && !isCurvedSegment) {
            // refit segment after recalibration
            TrkDriftCircleMath::DCSLFitter defaultFitter;
            TrkDriftCircleMath::Segment result(TrkDriftCircleMath::Line(0., 0., 0.), TrkDriftCircleMath::DCOnTrackVec());
            bool goodFit = defaultFitter.fit(result, line, segment.dcs(), hitSelector.selectHitsOnTrack(segment.dcs()));
            if (goodFit) {
                if (std::abs(xAOD::P4Helpers::deltaPhi(segment.line().phi(), result.line().phi())) > 0.01 ||
                    std::abs(segment.line().x0() - result.line().x0()) > 0.01 ||
                    std::abs(segment.line().y0() - result.line().y0()) > 0.01) {
                    // update local position and global
                    linephi = result.line().phi();
                    lpos[1] = result.line().position().x();
                    lpos[2] = result.line().position().y();
                    gpos = sInfo.amdbTrans * lpos;

                    // recreate  surface
                    surfaceTransform = Amg::Transform3D(sInfo.amdbTrans.rotation());
                    surfaceTransform.pretranslate(gpos);
                    surf = std::make_unique<Trk::PlaneSurface>(surfaceTransform, surfDim, surfDim);

                    // finally update global direction
                    gdir = updateDirection(linephi, *surf, roaddir2, isCurvedSegment);
                }
            }
        }

        // create local segment direction
        Trk::LocalDirection segLocDir;
        surf->globalToLocalDirection(gdir, segLocDir);
        if (segLocDir.angleYZ() == 0 && segLocDir.angleXZ() == 0) {
            ATH_MSG_DEBUG("invalid local direction");
            return nullptr;
        }

        // sanity checks
        const double diff_phi = xAOD::P4Helpers::deltaPhi(roaddir2.phi(), gdir.phi());
        const double diff_prec = xAOD::P4Helpers::deltaPhi(linephi, segLocDir.angleYZ());
        /// Use linearity of the sin at leading order to check that the angular differences are either 0 or PI
        if (std::min(std::abs(diff_phi), std::abs( std::abs(diff_phi) - M_PI)) > 1.e-3 || 
            std::min(std::abs(diff_prec), std::abs(std::abs(diff_prec) - M_PI)) > 1.e-3) {
            ATH_MSG_WARNING(" ALARM updated angles wrong: diff phi " << diff_phi << "  prec " << diff_prec << " phi rdir " << roaddir2.phi()
                                                                     << " gdir " << gdir.phi() << " lphi " << linephi << " seg "
                                                                     << segLocDir.angleYZ());
        }

        // associate Clusters to segment, uses spVecs to get clusters
        std::pair<std::pair<int, int>, bool> netaPhiHits =
            associateClustersToSegment(segment, chid, sInfo.globalTrans, sInfo.clusters, sInfo.phimin, sInfo.phimax, rioDistVec);

        if (rioDistVec.empty()){
            ATH_MSG_VERBOSE("No measurements were collected.");
            return nullptr;
            
        }
        /// Copy hits into vector
        auto meas_for_fit = [&rioDistVec] () {
            std::vector<const Trk::MeasurementBase*> out{};
            out.reserve(rioDistVec.size());
            std::sort(rioDistVec.begin(), rioDistVec.end(), SortByDistanceToSegment());
            for (const std::pair<double, std::unique_ptr<const Trk::MeasurementBase>>& ele : rioDistVec) out.push_back(ele.second.get());
            return out;
        };
        
            
        double dlocx{1000.}, dangleXZ{1000.}, qoverp{-99999.}, dqoverp{-99999.};
        bool hasMeasuredCoordinate = false;
        if (m_refitParameters && netaPhiHits.second) {
            ATH_MSG_DEBUG(" distance between first and last phi hit sufficient to perform 4D fit: phi  " << gdir.phi() << " theta "
                                                                                                         << gdir.theta());

            std::unique_ptr<Trk::Track> track{m_segmentFitter->fit(gpos, gdir, *surf, meas_for_fit())};

            if (track) {
                if (isCurvedSegment && track->perigeeParameters() && track->perigeeParameters()->covariance()) {
                    qoverp = track->perigeeParameters()->parameters()[Trk::qOverP];
                    dqoverp = Amg::error(*track->perigeeParameters()->covariance(), Trk::qOverP);
                }
                hasMeasuredCoordinate = true;
                // hack to update the second coordinate errors
                Amg::MatrixX updatedCov(5, 5);
                updatedCov.setZero();
                m_segmentFitter->updateSegmentParameters(*track, *surf, segLocPos, segLocDir, updatedCov);
                if (Amg::error(updatedCov, Trk::locX) > 0 && Amg::error(updatedCov, Trk::phi) > 0.) {
                    dlocx = Amg::error(updatedCov, Trk::locX);
                    dangleXZ = Amg::error(updatedCov, Trk::phi);  // hack (2): phi not always angleXZ
                } else {
                    ATH_MSG_WARNING(" Corrupt error matrix returned from fitter " << Amg::toString(updatedCov));
                }

                /// recalculate global direction and position
                surf->localToGlobal(segLocPos, gdir, gpos);
                surf->localToGlobalDirection(segLocDir, gdir);

                if (track->measurementsOnTrack() && rioDistVec.size() != track->measurementsOnTrack()->size()) {
                    if (track->measurementsOnTrack()->empty()) {
                        ATH_MSG_DEBUG("No measurements survived");
                        return nullptr;
                    }
                    ATH_MSG_DEBUG(" ROT vector size changed after fit, updating ");
                    garbage_collector = std::move(rioDistVec);
                    rioDistVec.reserve(track->measurementsOnTrack()->size());
                    const Trk::TrackParameters* firstPars = nullptr;
                    for (const Trk::TrackStateOnSurface* tsit : *track->trackStateOnSurfaces()) {
                        const Trk::TrackParameters* pars = tsit->trackParameters();
                        if (!pars) continue;
                        if (!firstPars) firstPars = pars;

                        // check whether state is a measurement, skip outliers if they are not MDT
                        const Trk::MeasurementBase* meas = tsit->measurementOnTrack();
                        if (!meas) continue;
                        if (tsit->type(Trk::TrackStateOnSurface::Outlier) && !dynamic_cast<const MdtDriftCircleOnTrack*>(meas)) continue;
                        double dist = (pars->position() - firstPars->position()).dot(firstPars->momentum().unit());
                        rioDistVec.emplace_back(dist, meas->uniqueClone());                        
                    }                    
                }
            } else {
                ATH_MSG_DEBUG(" refit of segment failed!! ");
                netaPhiHits.second = false;
            }
        }

        // optional update of phi position and direction, only performed if the segment was not refitted and there are phi
        // hits
        if (m_updatePhiUsingPhiHits && !netaPhiHits.second) {
            if (updateSegmentPhi(gpos, gdir, segLocPos, segLocDir, *surf, meas_for_fit(), sInfo.phimin, sInfo.phimax)) {
                surf->localToGlobal(segLocPos, gpos, gpos);
                surf->localToGlobalDirection(segLocDir, gdir);
                hasMeasuredCoordinate = true;
                dlocx = 100.;
                dangleXZ = 0.1;
            }
        }

        if (msgLvl(MSG::DEBUG)) {
            std::vector<const Trk::MeasurementBase*> debug_meas = meas_for_fit();
            ATH_MSG_DEBUG(" number of hits " << debug_meas.size() << " of which trigger " << netaPhiHits.first.first << " eta and "
                                             << netaPhiHits.first.second << " phi ");
            for (const Trk::MeasurementBase* mit : debug_meas) {
                const Trk::RIO_OnTrack* rot = dynamic_cast<const Trk::RIO_OnTrack*>(mit);
                if (rot) {
                    ATH_MSG_DEBUG(m_idHelperSvc->toString(rot->identify()));
                    const MdtDriftCircleOnTrack* mdt = dynamic_cast<const MdtDriftCircleOnTrack*>(rot);
                    if (mdt)
                        ATH_MSG_DEBUG(std::setprecision(4)
                                      << " radius " << std::setw(6) << mdt->driftRadius() << " time " << std::setw(6) << mdt->driftTime());
                    continue;
                }
                const CompetingMuonClustersOnTrack* crot = dynamic_cast<const CompetingMuonClustersOnTrack*>(mit);
                if (crot) {
                    ATH_MSG_DEBUG(m_idHelperSvc->toString(crot->rioOnTrack(0).identify())
                                  << " comp rot with hits " << crot->containedROTs().size());
                    continue;
                }
                ATH_MSG_WARNING("failed to dynamic_cast to ROT ");
            }
        }
        // recalculate holes
        std::vector<Identifier> holeVec = calculateHoles(ctx, chid, gpos, gdir, hasMeasuredCoordinate, deltaVec, outoftimeVec, rioDistVec);

        // currently not taking into account masked channels
        if (!outoftimeVec.empty()) holeVec.insert(holeVec.end(), outoftimeVec.begin(), outoftimeVec.end());
        MuonSegmentQuality* quality = new MuonSegmentQuality(segment.chi2(), segment.ndof(), holeVec);

        const TrkDriftCircleMath::DCSLFitter* dcslFitter = m_dcslFitProvider->getFitter();
        TrkDriftCircleMath::Segment result(TrkDriftCircleMath::Line(0., 0., 0.), TrkDriftCircleMath::DCOnTrackVec());
        if (dcslFitter && !segment.hasT0Shift() && m_outputFittedT0) {
            if (!dcslFitter->fit(result, segment.line(), segment.dcs(), hitSelector.selectHitsOnTrack(segment.dcs()))) {
                ATH_MSG_DEBUG(" T0 refit failed ");
            } else {              
                ATH_MSG_DEBUG(" Fitted T0 " << result.t0Shift()<<" is valid "<<result.hasT0Shift());                 
            }
        }
        bool hasFittedT0 = false;
        double fittedT0{0}, errorFittedT0{1.};
        if (m_outputFittedT0 && (segment.hasT0Shift() || (dcslFitter && result.hasT0Shift()))) {
            hasFittedT0 = true;
            if (segment.hasT0Shift()) {
                fittedT0 = segment.t0Shift();
                errorFittedT0 = segment.t0Error();
            } else if (dcslFitter && result.hasT0Shift()) {
                fittedT0 = result.t0Shift();
                errorFittedT0 = result.t0Error();
            } else {
                ATH_MSG_WARNING(" Failed to access fitted t0 ");
                hasFittedT0 = false;
            }
        }
        // create new segment
        std::unique_ptr<MuonSegment> msegment;
        if (isCurvedSegment) {  // curved segments
            if (qoverp == -99999.) {
                double charge = gpos.z() * std::tan(gdir.theta());
                charge = charge / std::abs(charge);
                // if the curved segment was not refit, then use a momentum estimate
                constexpr double BILALPHA(28.4366), BMLALPHA(62.8267), BMSALPHA(53.1259), BOLALPHA(29.7554);
                if (chIndex == MuonStationIndex::BIL) {
                    qoverp = (charge * segment.deltaAlpha()) / BILALPHA;
                    dqoverp = M_SQRT2 * segment.dtheta() / BILALPHA;
                } else if (chIndex == MuonStationIndex::BML) {
                    qoverp = (charge * segment.deltaAlpha()) / BMLALPHA;
                    dqoverp = M_SQRT2 * segment.dtheta() / BMLALPHA;
                } else if (chIndex == MuonStationIndex::BMS) {
                    qoverp = (charge * segment.deltaAlpha()) / BMSALPHA;
                    dqoverp = M_SQRT2 * segment.dtheta() / BMSALPHA;
                } else if (chIndex == MuonStationIndex::BOL) {
                    qoverp = (charge * segment.deltaAlpha()) / BOLALPHA;
                    dqoverp = M_SQRT2 * segment.dtheta() / BOLALPHA;
                }
            }
            Amg::MatrixX covMatrix(5, 5);
            covMatrix.setIdentity();
            covMatrix(0, 0) = dlocx * dlocx;
            covMatrix(1, 1) = segment.dy0() * segment.dy0();
            covMatrix(2, 2) = dangleXZ * dangleXZ;
            covMatrix(3, 3) = segment.dtheta() * segment.dtheta();
            covMatrix(4, 4) = dqoverp * dqoverp;

            std::vector<Trk::DefinedParameter> defPars;
            defPars.emplace_back(segLocPos[Trk::loc1], Trk::loc1);
            defPars.emplace_back(segLocPos[Trk::loc2], Trk::loc2);
            defPars.emplace_back(gdir.phi(), Trk::phi);
            defPars.emplace_back(gdir.theta(), Trk::theta);
            defPars.emplace_back(qoverp, Trk::qOverP);
            Trk::LocalParameters segLocPar(defPars);
            msegment = std::make_unique<MuonSegment>(segLocPar, covMatrix, surf.release(), createROTVec(rioDistVec).release(), quality,
                                                     Trk::Segment::DCMathSegmentMakerCurved);
        } else {  // straight segments
            // errors (for now no correlations)
            Amg::MatrixX covMatrix(4, 4);
            covMatrix.setIdentity();
            covMatrix(0, 0) = dlocx * dlocx;
            covMatrix(1, 1) = segment.dy0() * segment.dy0();
            covMatrix(2, 2) = dangleXZ * dangleXZ;
            covMatrix(3, 3) = segment.dtheta() * segment.dtheta();
            msegment = std::make_unique<MuonSegment>(segLocPos, segLocDir, covMatrix, surf.release(), createROTVec(rioDistVec).release(), quality,
                                                     Trk::Segment::DCMathSegmentMaker);
        }

        if (hasFittedT0) msegment->setT0Error(fittedT0, errorFittedT0);

        // check whether segment satisfies minimum quality criteria
        int segmentQuality = m_segmentSelectionTool->quality(*msegment);

        if (msgLvl(MSG::DEBUG)) {
            ATH_MSG_DEBUG(m_printer->print(*msegment) << " quality " << segmentQuality);
            if (segmentQuality < 0) ATH_MSG_DEBUG(" BAD segment ");
            if (hasFittedT0) ATH_MSG_DEBUG(" T0 " << fittedT0);
            if (isCurvedSegment) ATH_MSG_DEBUG(" Curved " << fittedT0);
        }
        if (segmentQuality < 0) { return nullptr; }
        return msegment;
    }

    void DCMathSegmentMaker::find(const std::vector<const Trk::RIO_OnTrack*>& rios, Trk::SegmentCollection* segColl) const {
        std::vector<const MdtDriftCircleOnTrack*> mdts;
        std::vector<const MuonClusterOnTrack*> clusters;

        for (const Trk::RIO_OnTrack* it : rios) {
            Identifier id = it->identify();
            if (m_idHelperSvc->isMdt(id)) {
                const MdtDriftCircleOnTrack* mdt = dynamic_cast<const MdtDriftCircleOnTrack*>(it);
                if (!mdt) { ATH_MSG_WARNING("failed dynamic_cast, not a MDT but hit has MDT id!!!"); }
                mdts.push_back(mdt);
            } else if (m_idHelperSvc->isTrigger(id)) {
                const MuonClusterOnTrack* clus = dynamic_cast<const MuonClusterOnTrack*>(it);
                if (!clus) { ATH_MSG_WARNING("failed dynamic_cast, not a cluster but hit has RPC/TGC id!!!"); }
                clusters.push_back(clus);
            }
        }
        find(mdts, clusters, segColl);
    }

    void DCMathSegmentMaker::find(const std::vector<const MdtDriftCircleOnTrack*>& mdts,
                                  const std::vector<const MuonClusterOnTrack*>& clusters, Trk::SegmentCollection* segColl) const {
        if (mdts.empty()) return;

        const MdtDriftCircleOnTrack* mdt = mdts.front();
        if (!mdt) return;

        bool hasPhiMeasurements = false;
        Amg::Vector3D gpos = mdt->globalPosition();
        Amg::Vector3D gdir = gpos.unit();
        find(gpos, gdir, mdts, clusters, hasPhiMeasurements, segColl);
    }

    void DCMathSegmentMaker::find(const std::vector<const Trk::RIO_OnTrack*>& rios1,
                                  const std::vector<const Trk::RIO_OnTrack*>& rios2) const {
        std::vector<const Trk::RIO_OnTrack*> rios = rios1;
        rios.insert(rios.end(), rios2.begin(), rios2.end());
        find(rios);
    }

    void DCMathSegmentMaker::find(const Trk::TrackRoad& road, const std::vector<std::vector<const MdtDriftCircleOnTrack*> >& mdts,
                                  const std::vector<std::vector<const MuonClusterOnTrack*> >& clusters, Trk::SegmentCollection* segColl,
                                  bool hasPhiMeasurements, double momentum) const {
        // copy all mdt hits into one vector
        std::vector<const MdtDriftCircleOnTrack*> all_mdts;
        for (const std::vector<const MdtDriftCircleOnTrack*>& circle_vec : mdts) { std::copy(circle_vec.begin(), circle_vec.end(), std::back_inserter(all_mdts)); }

        // copy all clusters into one vector
        std::vector<const MuonClusterOnTrack*> all_clus;
        for ( const std::vector<const MuonClusterOnTrack*>& clus_vec : clusters) { std::copy(clus_vec.begin(), clus_vec.end(), std::back_inserter(all_clus)); }

        const Amg::Vector3D& gpos = road.globalPosition();
        const Amg::Vector3D& gdir = road.globalDirection();
        find(gpos, gdir, all_mdts, all_clus, hasPhiMeasurements, segColl, momentum, road.deltaEta());
    }

    double DCMathSegmentMaker::errorScaleFactor(const Identifier& id, double curvature, bool hasPhiMeasurements) const {
        double scale = 1.;
        if (!m_curvedErrorScaling) return scale;

        if (!errorScalingRegion(id)) return scale;

        double scaleMax = 5.;
        if (m_curvedErrorScaling && curvature > 2) {
            scale = 1. + curvature / 10000;  // 3*(curvature/30000)
            if (scale > scaleMax) scale = scaleMax;
            ATH_MSG_DEBUG(" rescaled errors " << scale << " curvature " << curvature);
        }
        scale *= 2;

        // rescale errors is no phi measurement was found
        if (!hasPhiMeasurements) {
            double phiScale = 1.;
            // rescale errors
            int stRegion = m_idHelperSvc->mdtIdHelper().stationRegion(id);
            if (stRegion == 0)
                phiScale = 2.;  // inner
            else if (stRegion == 1)
                phiScale = 2.5;  // extended
            else if (stRegion == 2)
                phiScale = 2.5;  // middle
            else
                phiScale = 3.;  // outer
            scale = sqrt(scale * scale + phiScale * phiScale);
            ATH_MSG_DEBUG(" rescaled error for missing phi road " << scale);
        }

        return scale;
    }

    bool DCMathSegmentMaker::errorScalingRegion(const Identifier& id) const

    {
        // simple division of MuonSpectrometer in regions using barrel/endcap seperation plus
        // inner/middle/outer seperation

        bool isEndcap = m_idHelperSvc->isEndcap(id);

        if (isEndcap) {
            std::string stName = m_idHelperSvc->mdtIdHelper().stationNameString(m_idHelperSvc->mdtIdHelper().stationName(id));
            if (stName[1] == 'I') return true;

        } else {
            return true;
        }
        return false;
    }

    DCMathSegmentMaker::ClusterVecPair DCMathSegmentMaker::create1DClusters(const std::vector<const MuonClusterOnTrack*>& clusters) const {
        // if empty return
        if (clusters.empty()) return ClusterVecPair();
        // some useful typedefs...
        typedef std::vector<const MuonClusterOnTrack*>::const_iterator ROTCit;

        // create a vector to hold the clusters
        ClusterVec clVec;
        ClusterVec phiVec;
        clVec.reserve(clusters.size());

        ROTCit cit = clusters.begin();
        ROTCit cit_end = clusters.end();
        for (; cit != cit_end; ++cit) {
            Identifier id = (*cit)->identify();
            Identifier gasGapId = m_idHelperSvc->gasGapId(id);
            int measuresPhi = m_idHelperSvc->measuresPhi(id);

            if (measuresPhi) {
                phiVec.push_back(createSpacePoint(gasGapId, nullptr, *cit));
                if (phiVec.back().corrupt()) phiVec.pop_back();
            } else {
                clVec.push_back(createSpacePoint(gasGapId, *cit, nullptr));
                if (clVec.back().corrupt()) clVec.pop_back();
            }
        }

        return ClusterVecPair(clVec, phiVec);
    }

    DCMathSegmentMaker::ClusterVecPair DCMathSegmentMaker::create2DClusters(const std::vector<const MuonClusterOnTrack*>& clusters) const {
        // if empty return
        if (clusters.empty()) return ClusterVecPair();

        ChIdHitMap gasGapHitMap;

        ROTCit cit = clusters.begin();
        ROTCit cit_end = clusters.end();
        for (; cit != cit_end; ++cit) {
            Identifier id = (*cit)->identify();
            ATH_MSG_VERBOSE(" new trigger hit " << m_idHelperSvc->toString(id));

            Identifier chId = m_idHelperSvc->chamberId(id);
            Identifier gasGapId = m_idHelperSvc->gasGapId(id);
            int measuresPhi = m_idHelperSvc->measuresPhi(id);

            // eta hits first than phi hits
            if (!measuresPhi)
                gasGapHitMap[chId][gasGapId].first.push_back(*cit);
            else
                gasGapHitMap[chId][gasGapId].second.push_back(*cit);
        }

        return createSpacePoints(gasGapHitMap);
    }

    void 
    DCMathSegmentMaker::handleChamber(DCMathSegmentMaker::IdHitMap& gasGapHitMap) const {
        ATH_MSG_DEBUG("  new chamber " << m_idHelperSvc->toString(gasGapHitMap.begin()->first));
        std::list<const Trk::PrepRawData*> prds;

        IdHitMapIt git = gasGapHitMap.begin();
        IdHitMapIt git_end = gasGapHitMap.end();
        for (; git != git_end; ++git) {
            ROTCit hit = git->second.first.begin();
            ROTCit hit_end = git->second.first.end();
            for (; hit != hit_end; ++hit) {
                double lpos = (*hit)->localParameters()[Trk::locX];
                ATH_MSG_DEBUG("     new eta hit  " << lpos << "  " << m_idHelperSvc->toString((*hit)->identify()));
                prds.push_back((*hit)->prepRawData());
            }
        }
        std::unique_ptr<const CompetingMuonClustersOnTrack> rotEta;
        if (prds.size() > 1) {
            rotEta = m_compClusterCreator->createBroadCluster(prds, 0.);

            ATH_MSG_DEBUG("     --> new competing eta rot  " << rotEta->localParameters()[Trk::locX] << " "
                                                             << Amg::error(rotEta->localCovariance(), Trk::locX));
        }
        prds.clear();
        git = gasGapHitMap.begin();
        git_end = gasGapHitMap.end();
        for (; git != git_end; ++git) {
            ROTCit hit = git->second.second.begin();
            ROTCit hit_end = git->second.second.end();
            for (; hit != hit_end; ++hit) {
                double lpos = (*hit)->localParameters()[Trk::locX];

                ATH_MSG_DEBUG("     new phi hit  " << lpos << "  " << m_idHelperSvc->toString((*hit)->identify()));
                prds.push_back((*hit)->prepRawData());
            }
        }
        std::unique_ptr<const CompetingMuonClustersOnTrack> rotPhi;
        if (prds.size() > 1) {
            rotPhi = m_compClusterCreator->createBroadCluster(prds, 0.);

            ATH_MSG_DEBUG("     --> new competing phi rot  " << rotPhi->localParameters()[Trk::locX] << " "
                                                             << Amg::error(rotPhi->localCovariance(), Trk::locX));
        }
    }

    DCMathSegmentMaker::ClusterVecPair DCMathSegmentMaker::createSpacePoints(DCMathSegmentMaker::ChIdHitMap& chIdHitMap) const {
        // vector to store output
        ClusterVecPair clusterVecs;
        clusterVecs.first.reserve(20);
        clusterVecs.second.reserve(20);

        // loop over chambers
        ChIdHitMapIt cit = chIdHitMap.begin();
        ChIdHitMapIt cit_end = chIdHitMap.end();
        for (; cit != cit_end; ++cit) {
            // create clusters per chamber and copy them in to result vector
            ClusterVecPair cls = createSpacePoints(cit->second);
            std::copy(cls.first.begin(), cls.first.end(), std::back_inserter(clusterVecs.first));
            std::copy(cls.second.begin(), cls.second.end(), std::back_inserter(clusterVecs.second));

            // handleChamber( cit->second );
        }

        return clusterVecs;
    }

    DCMathSegmentMaker::ClusterVecPair DCMathSegmentMaker::createSpacePoints(DCMathSegmentMaker::IdHitMap& gasGapHitMap) const {
        ClusterVec spacePoints;
        ClusterVec phiVec;
        bool isEndcap = m_idHelperSvc->isEndcap((*(gasGapHitMap.begin())).first);

        ATH_MSG_VERBOSE(" creating Space points for " << gasGapHitMap.size() << " gas gaps ");

        IdHitMapIt git = gasGapHitMap.begin();
        IdHitMapIt git_end = gasGapHitMap.end();
        for (; git != git_end; ++git) {
            // flag whether phi hits are matched with a eta hit
            std::vector<int> flagPhihit(git->second.second.size(), 0);

            // store Identifier of previous hit to remove duplicates
            Identifier prevEtaId;

            ATH_MSG_VERBOSE(" New gasgap " << m_idHelperSvc->toString(git->first) << " neta " << git->second.first.size() << " nphi "
                                           << git->second.second.size());

            ROTCit etaIt = git->second.first.begin();
            ROTCit etaIt_end = git->second.first.end();
            for (; etaIt != etaIt_end; ++etaIt) {
                // check whether we are not dealing with a duplicate hit
                if ((*etaIt)->identify() == prevEtaId) continue;
                prevEtaId = (*etaIt)->identify();

                ATH_MSG_VERBOSE(" Eta hit " << m_idHelperSvc->toString((*etaIt)->identify()));

                if (isEndcap) {
                    // check whether match with phi hits was found
                    bool foundSP = false;
                    Identifier prevPhiId;
                    ROTCit phiIt = git->second.second.begin();
                    ROTCit phiIt_begin = git->second.second.begin();
                    ROTCit phiIt_end = git->second.second.end();
                    for (; phiIt != phiIt_end; ++phiIt) {
                        // check for duplicate phi hits
                        if ((*phiIt)->identify() == prevPhiId) continue;
                        prevPhiId = (*phiIt)->identify();

                        ATH_MSG_VERBOSE(" Phi hit " << m_idHelperSvc->toString((*phiIt)->identify()));

                        Cluster2D sp = createTgcSpacePoint(git->first, *etaIt, *phiIt);
                        if (sp.corrupt()) continue;
                        spacePoints.push_back(sp);

                        // mark as used
                        foundSP = true;
                        flagPhihit[std::distance(phiIt_begin, phiIt)] = 1;
                    }

                    // add single eta hit if not matching phi hit was found
                    if (!foundSP) {
                        Cluster2D sp = createSpacePoint(git->first, *etaIt, nullptr);
                        if (sp.corrupt()) continue;
                        spacePoints.push_back(sp);
                    }
                } else {
                    Cluster2D sp = createRpcSpacePoint(git->first, *etaIt, git->second.second);
                    if (sp.corrupt()) continue;
                    // flag all phi hits, not very elegant, but works
                    flagPhihit = std::vector<int>(git->second.second.size(), 1);
                    spacePoints.push_back(sp);
                }
            }
            if (isEndcap) {
                // loop over flag vector and add unmatched phi hits to phiVec;
                Identifier prevPhiId;
                for (unsigned int i = 0; i < flagPhihit.size(); ++i) {
                    if (flagPhihit[i]) continue;

                    // check for duplicate phi hits
                    if (git->second.second[i]->identify() == prevPhiId) continue;
                    prevPhiId = git->second.second[i]->identify();

                    Cluster2D sp = createTgcSpacePoint(git->first, nullptr, git->second.second[i]);
                    if (sp.corrupt()) continue;
                    phiVec.push_back(sp);
                }
            } else if (git->second.first.empty() && !git->second.second.empty()) {
                // if there were no eta hits create one phi spacePoint of all phi hits in gasgap
                Cluster2D sp = createRpcSpacePoint(git->first, nullptr, git->second.second);
                if (sp.corrupt()) continue;
                phiVec.push_back(sp);
            }
        }

        ATH_MSG_VERBOSE(" Creating space points, number of gas-gaps " << gasGapHitMap.size() << "  space points " << spacePoints.size());

        return ClusterVecPair(spacePoints, phiVec);
    }

    DCMathSegmentMaker::Cluster2D DCMathSegmentMaker::createSpacePoint(const Identifier& gasGapId, const MuonClusterOnTrack* etaHit,
                                                                       const MuonClusterOnTrack* phiHit) const {
        bool isEndcap = m_idHelperSvc->isEndcap(gasGapId);
        double error{1.}, lpx{0.}, lpy{0.};
        // case one hit missing. Take position and error of the available hit
        if (!etaHit) {
            if (!phiHit) {
                ATH_MSG_WARNING("Both eta and phi hits missing");
                lpx = 0;
                error = 0;
            } else {
                lpx = phiHit->localParameters()[Trk::locX];
                error = Amg::error(phiHit->localCovariance(), Trk::locX);
            }
        } else if (!phiHit) {
            lpx = etaHit->localParameters()[Trk::locX];
            error = Amg::error(etaHit->localCovariance(), Trk::locX);
        } else if (etaHit && phiHit) {
            if (isEndcap) {
                return createTgcSpacePoint(gasGapId, etaHit, phiHit);
            } else {
                std::vector<const MuonClusterOnTrack*> phiVec;
                phiVec.push_back(phiHit);
                return createRpcSpacePoint(gasGapId, etaHit, phiVec);
            }
        }
        Identifier detElId = m_idHelperSvc->detElId(gasGapId);
        if (std::abs(error) < 0.001) {
            ATH_MSG_WARNING(" Unphysical error assigned for gasgap " << m_idHelperSvc->toString(gasGapId));
            error = 0.;
        }
        return Cluster2D(detElId, gasGapId, Amg::Vector2D(lpx, lpy), error, etaHit, phiHit);
    }

    DCMathSegmentMaker::Cluster2D DCMathSegmentMaker::createTgcSpacePoint(const Identifier& gasGapId, const MuonClusterOnTrack* etaHit,
                                                                          const MuonClusterOnTrack* phiHit) const {
        double error{1.}, lpx{0.}, lpy{0.};
        Identifier detElId = m_idHelperSvc->detElId(gasGapId);
        // case one hit missing. Take position and error of the available hit
        if (!etaHit) {
            lpx = phiHit->localParameters()[Trk::locX];
            error = Amg::error(phiHit->localCovariance(), Trk::locX);
        } else if (!phiHit) {
            lpx = etaHit->localParameters()[Trk::locX];
            error = Amg::error(etaHit->localCovariance(), Trk::locX);
        } else if (etaHit && phiHit) {
            // get orientation angle of strip to rotate back from local frame to strip
            // copy code from ROT creator
            int stripNo = m_idHelperSvc->tgcIdHelper().channel(phiHit->identify());
            int gasGap = m_idHelperSvc->tgcIdHelper().gasGap(phiHit->identify());

            const MuonGM::TgcReadoutElement* detEl = dynamic_cast<const MuonGM::TgcReadoutElement*>(etaHit->detectorElement());
            if (!detEl) {
                ATH_MSG_WARNING("dynamic cast error for "<<m_idHelperSvc->toString(etaHit->identify())<<". Expected TGCs. Returning");
                return Cluster2D(detElId, gasGapId, Amg::Vector2D(lpx, lpy), error, etaHit, phiHit);
            }
            // calculate local position of endpoint of strip
            Amg::Vector3D lEtapos = detEl->localChannelPos(etaHit->identify());
            double localEtaY = detEl->stripCtrX(gasGap, stripNo, lEtapos.z());
            if (0 < detEl->getStationEta()) { localEtaY *= -1.; }
            Amg::Vector3D lSppos = lEtapos;
            lSppos[1] = localEtaY;

            // transform to global
            const Amg::Transform3D tgcTrans = detEl->absTransform();
            Amg::Vector3D gposSp = tgcTrans * lSppos;
            lpx = etaHit->localParameters()[Trk::locX];
            error = Amg::error(etaHit->localCovariance(), Trk::locX);
            if (error == 0.) {
                ATH_MSG_WARNING(" Unphysical error assigned for " << m_idHelperSvc->toString(etaHit->identify()));
                if (etaHit->prepRawData())
                    ATH_MSG_WARNING(" PRD error " << Amg::error(etaHit->prepRawData()->localCovariance(), Trk::locX));
            }
            Amg::Vector2D lspPos(0., 0.);
            if (etaHit->associatedSurface().globalToLocal(gposSp, gposSp, lspPos)) {
                lpy = lspPos[Trk::locY];
            } else {
                ATH_MSG_WARNING(" globalToLocal failed ");
            }

            ATH_MSG_DEBUG(" TGC space point: error " << error << " stripWith " << error * M_SQRT2 << std::endl
                                                     << "   " << m_idHelperSvc->toString(etaHit->identify()) << std::endl
                                                     << "   " << m_idHelperSvc->toString(phiHit->identify()));
        }
        if (std::abs(error) < 0.001) {
            ATH_MSG_WARNING(" Unphysical error assigned for gasgap " << m_idHelperSvc->toString(gasGapId));
            error = 1.;
        }
        return Cluster2D(detElId, gasGapId, Amg::Vector2D(lpx, lpy), error, etaHit, phiHit);
    }

    DCMathSegmentMaker::Cluster2D DCMathSegmentMaker::createRpcSpacePoint(const Identifier& gasGapId, const MuonClusterOnTrack* etaHit,
                                                                          const std::vector<const MuonClusterOnTrack*>& phiHits) const {
        // create vector to store phi hits after removal of duplicate hits
        std::vector<const MuonClusterOnTrack*> cleanPhihits;
        cleanPhihits.reserve(phiHits.size());

        double error{1.}, lpx{0.}, lpy{0.};
        // case one hit missing. Take position and error of the available hit
        if (!etaHit) {
            lpx = phiHits.front()->localParameters()[Trk::locX];
            error = Amg::error(phiHits.front()->localCovariance(), Trk::locX);
            // loop over phi hits, remove duplicate phi hits
            Identifier prevId;
            for (const MuonClusterOnTrack* clus : phiHits) {
                // remove duplicate phi hits
                if (clus->identify() == prevId) continue;
                prevId = clus->identify();
                cleanPhihits.push_back(clus);
            }
        } else if (phiHits.empty()) {
            lpx = etaHit->localParameters()[Trk::locX];
            error = Amg::error(etaHit->localCovariance(), Trk::locX);
        } else if (etaHit && !phiHits.empty()) {
            lpx = etaHit->localParameters()[Trk::locX];
            error = Amg::error(etaHit->localCovariance(), Trk::locX);

            ATH_MSG_DEBUG(" RPC space point: error " << error << " stripWith " << error * M_SQRT2 << std::endl
                                                     << "   " << m_idHelperSvc->toString(etaHit->identify()));

            double minPos{1e9}, maxPos{-1e9};
            Identifier prevId;

            // loop over phi hits, calculate average position + cluster width, remove duplicate phi hits
            std::vector<const MuonClusterOnTrack*>::const_iterator pit = phiHits.begin();
            std::vector<const MuonClusterOnTrack*>::const_iterator pit_end = phiHits.end();
            for (; pit != pit_end; ++pit) {
                // remove duplicate phi hits
                if ((*pit)->identify() == prevId) continue;
                prevId = (*pit)->identify();

                // calculate phi hit position in local eta hit reference frame
                Amg::Vector2D phiLocPos(0., 0.);
                if (etaHit->associatedSurface().globalToLocal((*pit)->globalPosition(), (*pit)->globalPosition(), phiLocPos)) {
                    lpy = phiLocPos[Trk::locY];
                    if (lpy < minPos) minPos = lpy;
                    if (lpy > maxPos) maxPos = lpy;
                    ATH_MSG_DEBUG("    " << m_idHelperSvc->toString((*pit)->identify()));
                    cleanPhihits.push_back(*pit);
                }
            }
            if (cleanPhihits.size() > 1)
                ATH_MSG_DEBUG("  multiple phi hits: nhits " << cleanPhihits.size() << " cl width " << maxPos - minPos);
        } else {
            ATH_MSG_DEBUG(" ARRRGGG got two empty pointers!!! ");
        }
        Identifier detElId = m_idHelperSvc->detElId(gasGapId);
        if (std::abs(error) < 0.001) {
            ATH_MSG_WARNING(" Unphysical error assigned for gasgap " << m_idHelperSvc->toString(gasGapId));
            error = 1.;
        }
        return Cluster2D(detElId, gasGapId, Amg::Vector2D(lpx, lpy), error, etaHit, !cleanPhihits.empty() ? cleanPhihits : phiHits);
    }

    TrkDriftCircleMath::CLVec DCMathSegmentMaker::createClusterVec(const Identifier& chid, ClusterVec& spVec,
                                                                   const Amg::Transform3D& gToStation) const {
        TrkDriftCircleMath::CLVec cls;

        int chPhi = m_idHelperSvc->mdtIdHelper().stationPhi(chid);

        // loop over clusters
        ClusterIt cit = spVec.begin();
        ClusterIt cit_end = spVec.end();
        int index = 0;
        cls.reserve(spVec.size());
        for (; cit != cit_end; ++cit) {
            const MuonClusterOnTrack* meas = cit->etaHit ? cit->etaHit : cit->phiHit;
            // construct cluster id
            Identifier id = meas->identify();
            int measuresPhi = m_idHelperSvc->measuresPhi(id);
            int eta = m_idHelperSvc->stationEta(id);
            int phi = m_idHelperSvc->stationPhi(id);
            int isTgc = m_idHelperSvc->isTgc(id);
            int name = isTgc ? m_idHelperSvc->tgcIdHelper().stationName(id) : m_idHelperSvc->rpcIdHelper().stationName(id);
            if (!isTgc) {
                if (chPhi != phi) {
                    ATH_MSG_VERBOSE(" Discarding cluster, wrong station phi " << m_idHelperSvc->toString(id));
                    ++index;
                    continue;
                }
            }
            TrkDriftCircleMath::ClusterId clid(name, eta, phi, isTgc, measuresPhi);

            // calculate local cluster position
            Amg::Vector3D locPos = gToStation * cit->globalPos;
            TrkDriftCircleMath::LocVec2D lp(locPos.y(), locPos.z());

            if (std::abs(lp.y()) > m_maxAssociateClusterDistance) {
                ATH_MSG_VERBOSE(" Discarding cluster with large distance from chamber " << m_idHelperSvc->toString(id));
                ++index;
                continue;
            }
            ATH_MSG_VERBOSE(" " << m_idHelperSvc->toString(id) << "  clid: " << clid.id() << " central phi "
                                << meas->detectorElement()->center().phi() << " index " << index);
            double err = cit->error;
            TrkDriftCircleMath::Cluster cl(lp, err, clid, index);
            cls.push_back(cl);

            ++index;
        }
        return cls;
    }

    TrkDriftCircleMath::DCVec DCMathSegmentMaker::createDCVec(const std::vector<const MdtDriftCircleOnTrack*>& mdts, double errorScale,
                                                              std::set<Identifier>& chamberSet, double& phimin, double& phimax,
                                                              TrkDriftCircleMath::DCStatistics& dcStatistics, const Amg::Transform3D& gToStation,
                                                              const Amg::Transform3D& amdbToGlobal) const {
        TrkDriftCircleMath::DCVec dcs;
        dcs.reserve(mdts.size());
        /* ********  Mdt hits  ******** */

        //     typedef std::map< Identifier, TrkDriftCircleMath::MdtChamberGeometry > GeoMap;
        //     typedef GeoMap::iterator GeoIt;
        //     GeoMap geoMap;

        bool firstMdt = true;

        unsigned index = 0;
        for (const MdtDriftCircleOnTrack* rot : mdts) {
            /// Need to chech the constructor of Drift Circle. Index does not seem to be set properly
            if (!rot) {
                ATH_MSG_WARNING(" rot not a MdtDriftCircleOnTrack ");
                ++index;
                continue;
            }
            const MuonGM::MdtReadoutElement* detEl = rot->prepRawData()->detectorElement();

            if (!detEl) {
                ATH_MSG_WARNING(" aborting not detEl found ");
                return TrkDriftCircleMath::DCVec();
            }

            Identifier id = rot->identify();
            Identifier elId = m_idHelperSvc->mdtIdHelper().elementID(id);

            // calculate local AMDB position
            Amg::Vector3D locPos = gToStation * rot->prepRawData()->globalPosition();
            TrkDriftCircleMath::LocVec2D lpos(locPos.y(), locPos.z());

            double r = rot->localParameters()[Trk::locR];
            double dr = Amg::error(rot->localCovariance(), Trk::locR) * errorScale;

            // create identifier
            TrkDriftCircleMath::MdtId mdtid(m_idHelperSvc->mdtIdHelper().isBarrel(id), m_idHelperSvc->mdtIdHelper().multilayer(id) - 1,
                                            m_idHelperSvc->mdtIdHelper().tubeLayer(id) - 1, m_idHelperSvc->mdtIdHelper().tube(id) - 1);

            //
            double preciseError = dr;
            if (m_usePreciseError) { preciseError = m_preciseErrorScale * (0.23 * std::exp(-std::abs(r) / 6.06) + 0.0362); }
            // create new DriftCircle
            TrkDriftCircleMath::DriftCircle dc(lpos, r, dr, preciseError, TrkDriftCircleMath::DriftCircle::InTime, std::move(mdtid), index, rot);

            TubeEnds tubeEnds = localTubeEnds(*rot, gToStation, amdbToGlobal);
            if (firstMdt) {
                phimin = tubeEnds.phimin;
                phimax = tubeEnds.phimax;
                firstMdt = false;
            } else {
                updatePhiRanges(tubeEnds.phimin, tubeEnds.phimax, phimin, phimax);
            }

            if (msgLvl(MSG::VERBOSE)) {
                ATH_MSG_VERBOSE(" new MDT hit " << m_idHelperSvc->toString(id) << " x " << lpos.x() << " y " << lpos.y() << " time "
                                                << rot->driftTime() << " r " << r << " dr " << dr << " phi range " << tubeEnds.phimin << " "
                                                << tubeEnds.phimax);
                if (m_usePreciseError) ATH_MSG_VERBOSE(" dr(2) " << preciseError);
            }
            dcs.push_back(std::move(dc));

            chamberSet.insert(elId);

            ++dcStatistics[detEl];

            ++index;
        }

        return dcs;
    }

    TrkDriftCircleMath::MdtChamberGeometry DCMathSegmentMaker::createChamberGeometry(const Identifier& chid,
                                                                                     const Amg::Transform3D& gToStation) const {
        /* calculate chamber geometry
           it takes as input:
           distance between the first and second tube in the chamber within a layer along the tube layer (tube distance)
           distance between the first tube in the first layer and the first tube in the second layer along the tube layer
           (tube stagering) distance between the first and second layer perpendicular to the tube layers (layer distance)
           position of the first hit in ml 0 and ml 1 (2D in plane)
           total number of multilayers
           total number of layers
           total number of tubes per layer for each multilayer
           an identifier uniquely identifying the chamber
        */

        Amg::Vector3D firstTubeMl0(0., 0., 0.);
        Amg::Vector3D firstTubeMl1(0., 0., 0.);

        // get id
        int eta = m_idHelperSvc->mdtIdHelper().stationEta(chid);
        int phi = m_idHelperSvc->mdtIdHelper().stationPhi(chid);
        int name = m_idHelperSvc->mdtIdHelper().stationName(chid);
        
        SG::ReadCondHandle<MuonGM::MuonDetectorManager> DetectorManagerHandle{m_DetectorManagerKey};
        const MuonGM::MuonDetectorManager* MuonDetMgr{*DetectorManagerHandle};
        if (!MuonDetMgr) { 
          ATH_MSG_ERROR("Null pointer to the read MuonDetectorManager conditions object"); 
          return TrkDriftCircleMath::MdtChamberGeometry();
        }

        // get detEL for first ml (always there)
        const MuonGM::MdtReadoutElement* detEl1 =
            MuonDetMgr->getMdtReadoutElement(m_idHelperSvc->mdtIdHelper().channelID(name, eta, phi, 1, 1, 1));
        const MuonGM::MdtReadoutElement* detEl2 = nullptr;
        int ntube2 = 0;
        // number of multilayers in chamber
        int nml = detEl1->nMDTinStation();

        // treament of chambers with two ml
        if (nml == 2) {
            Identifier firstIdml1 = m_idHelperSvc->mdtIdHelper().channelID(name, eta, phi, 2, 1, 1);
            detEl2 = MuonDetMgr->getMdtReadoutElement(firstIdml1);
            firstTubeMl1 = gToStation * (detEl2->surface(firstIdml1).center());
            ntube2 = detEl2->getNtubesperlayer();
        }

        // number of layers and tubes
        int nlay = detEl1->getNLayers();
        int ntube1 = detEl1->getNtubesperlayer();

        // position first tube in ml 0
        Identifier firstIdml0 = m_idHelperSvc->mdtIdHelper().channelID(name, eta, phi, 1, 1, 1);
        firstTubeMl0 = gToStation * (detEl1->surface(firstIdml0).center());

        // position second tube in ml 0
        Identifier secondIdml0 = m_idHelperSvc->mdtIdHelper().channelID(name, eta, phi, 1, 1, 2);
        Amg::Vector3D secondTubeMl0 = gToStation * (detEl1->surface(secondIdml0).center());

        TrkDriftCircleMath::LocVec2D firstTube0(firstTubeMl0.y(), firstTubeMl0.z());
        TrkDriftCircleMath::LocVec2D firstTube1(firstTubeMl1.y(), firstTubeMl1.z());

        // position first tube ml 0 and 1
        Identifier firstIdml0lay1 = m_idHelperSvc->mdtIdHelper().channelID(name, eta, phi, 1, 2, 1);
        Amg::Vector3D firstTubeMl0lay1 = gToStation * (detEl1->surface(firstIdml0lay1).center());

        double tubeDist = (secondTubeMl0 - firstTubeMl0).y();      // distance between tube in a given layer
        double tubeStage = (firstTubeMl0lay1 - firstTubeMl0).y();  // tube stagering distance
        double layDist = (firstTubeMl0lay1 - firstTubeMl0).z();    // distance between layers

        TrkDriftCircleMath::MdtChamberGeometry mdtgeo(chid, m_idHelperSvc.get(), nml, nlay, ntube1, ntube2, firstTube0, firstTube1, tubeDist, tubeStage,
                                                      layDist, detEl1->surface().center().theta());

        if (msgLvl(MSG::VERBOSE)) mdtgeo.print(msgStream());

        return mdtgeo;
    }

   void DCMathSegmentMaker::associateMDTsToSegment(
        const Amg::Vector3D& gdir, TrkDriftCircleMath::Segment& segment, const std::vector<const MdtDriftCircleOnTrack*>& mdts,
        TrkDriftCircleMath::MdtMultiChamberGeometry* multiGeo, const Amg::Transform3D& gToStation, const Amg::Transform3D& amdbToGlobal,
        std::set<Identifier>& deltaVec, std::set<Identifier>& outoftimeVec,
        std::vector<std::pair<double,  std::unique_ptr<const Trk::MeasurementBase>> >& rioDistVec) const {
        // clear result vectors
    
        // convert segment parameters + x position from road
        const TrkDriftCircleMath::Line& line = segment.line();
        TrkDriftCircleMath::TransformToLine toLineml1(line);
        TrkDriftCircleMath::TransformToLine toLineml2(line);
        if (segment.hasCurvatureParameters()) {
            // ml2 segment direction
            double ml2phi = line.phi() - segment.deltaAlpha();
            TrkDriftCircleMath::LocVec2D ml2dir(std::cos(ml2phi), std::sin(ml2phi));
            // ml2 segment position
            const TrkDriftCircleMath::LocVec2D ml1LocPos = multiGeo->tubePosition(0, multiGeo->nlay(), 0);
            const TrkDriftCircleMath::LocVec2D ml2LocPos = multiGeo->tubePosition(1, 1, 0);
            double chamberMidPtY = (ml1LocPos.y() + ml2LocPos.y()) / 2.0;
            TrkDriftCircleMath::LocVec2D ml2pos(segment.deltab(), chamberMidPtY);
            // construct the new ml2 segment line & transform
            const TrkDriftCircleMath::Line ml2line(ml2pos, ml2dir);
            TrkDriftCircleMath::TransformToLine tmptoLine(ml2line);
            // set the ml2 line
            toLineml2 = tmptoLine;
        }

        for (TrkDriftCircleMath::DCOnTrack& dcit : segment.dcs()) {
            if (dcit.state() == TrkDriftCircleMath::DCOnTrack::Delta) { deltaVec.insert(mdts[dcit.index()]->identify()); }

            if (dcit.state() == TrkDriftCircleMath::DCOnTrack::OutOfTime) { outoftimeVec.insert(mdts[dcit.index()]->identify()); }

            if (dcit.state() != TrkDriftCircleMath::DCOnTrack::OnTrack) continue;

            const MdtDriftCircleOnTrack* riodc = dynamic_cast<const MdtDriftCircleOnTrack*>(mdts[dcit.index()]);
            if (!riodc) continue;

            // choose which line to use (ml1 or ml2)
            TrkDriftCircleMath::TransformToLine toLine = toLineml1;
            if (m_idHelperSvc->mdtIdHelper().multilayer(riodc->identify()) == 2) toLine = toLineml2;
            // calculate position of hit in line frame
            TrkDriftCircleMath::LocVec2D pointOnHit = toLine.toLine(dcit.position());

            // calculate position of hit on line in line frame
            TrkDriftCircleMath::LocVec2D pointOnLine(pointOnHit.x(), 0.);

            // transform back to local AMDB coordinates
            TrkDriftCircleMath::LocVec2D pointOnLineAMDB = toLine.toLocal(pointOnLine);

            // get position along wire from ROT
            Amg::Vector3D posAlong = gToStation * riodc->globalPosition();

            // set yz components
            posAlong[1] = pointOnLineAMDB.x();
            posAlong[2] = pointOnLineAMDB.y();

            // back to global
            Amg::Vector3D mdtGP = amdbToGlobal * posAlong;

            const Trk::StraightLineSurface* surf = dynamic_cast<const Trk::StraightLineSurface*>(&riodc->associatedSurface());
            if (!surf) {
                ATH_MSG_WARNING(" dynamic cast to StraightLineSurface failed for mdt!!! ");
                continue;
            }

            // calculate Amg::Vector2D using surf to obtain sign
            Amg::Vector2D locPos(0., 0.);
            if (!surf->globalToLocal(mdtGP, gdir, locPos)) ATH_MSG_WARNING(" globalToLocal failed ");

            // calculate side
            Trk::DriftCircleSide side = locPos[Trk::driftRadius] < 0 ? Trk::LEFT : Trk::RIGHT;

            std::unique_ptr<MdtDriftCircleOnTrack> nonconstDC;
            bool hasT0 = segment.hasT0Shift();
            if (!hasT0) {
                // ATH_MSG_VERBOSE(" recalibrate MDT hit");
                nonconstDC.reset(m_mdtCreator->createRIO_OnTrack(*riodc->prepRawData(), mdtGP, &gdir));
            } else {
                ATH_MSG_VERBOSE(" recalibrate MDT hit with shift " << segment.t0Shift());
                nonconstDC.reset(m_mdtCreatorT0->createRIO_OnTrack(*riodc->prepRawData(), mdtGP, &gdir, segment.t0Shift()));
            }

            if (!nonconstDC) {
                dcit.state(TrkDriftCircleMath::DCOnTrack::OutOfTime);
                continue;
            }

            // update the drift radius after recalibration, keep error
            TrkDriftCircleMath::DriftCircle new_dc(dcit.position(), std::abs(nonconstDC->driftRadius()), dcit.dr(), dcit.drPrecise(),
                                                   dcit.driftState(), dcit.id(),
                                                   dcit.index(), nonconstDC.get());
            TrkDriftCircleMath::DCOnTrack new_dc_on_track(std::move(new_dc), dcit.residual(), dcit.errorTrack());
            dcit = std::move(new_dc_on_track);

            if (hasT0) {
                if (msgLvl(MSG::VERBOSE)) {
                    double shift = riodc->driftTime() - nonconstDC->driftTime();
                    ATH_MSG_VERBOSE(" t0 shift " << segment.t0Shift() << " from hit " << shift << " recal " << nonconstDC->driftRadius()
                                                 << " t " << nonconstDC->driftTime() << "  from fit " << dcit.r() << " old "
                                                 << riodc->driftRadius() << " t " << riodc->driftTime());
                    if (std::abs(std::abs(nonconstDC->driftRadius()) - std::abs(dcit.r())) > 0.1 && nonconstDC->driftRadius() < 19. &&
                        nonconstDC->driftRadius() > 1.) {
                        ATH_MSG_WARNING("Detected invalid recalibration after T0 shift");
                    }
                }
            }
            m_mdtCreator->updateSign(*nonconstDC, side);
            double dist = pointOnHit.x();
            rioDistVec.push_back(std::make_pair(dist, std::move(nonconstDC)));
        }       
    }

    std::string DCMathSegmentMaker::printSP(std::pair<double, double> resPull, const DCMathSegmentMaker::Cluster2D& spacePoint) const {
        std::ostringstream sout;
        Identifier id;
        if (spacePoint.etaHit != nullptr && spacePoint.phiHit != nullptr) {
            sout << " space point: ";
            id = spacePoint.etaHit->identify();
        } else if (spacePoint.etaHit != nullptr && spacePoint.phiHit == nullptr) {
            sout << " eta hit:     ";
            id = spacePoint.etaHit->identify();
        } else if (spacePoint.etaHit == nullptr && spacePoint.phiHit != nullptr) {
            sout << " phi hit:     ";
            id = spacePoint.phiHit->identify();
        } else
            sout << " UNKNOWN TYPE ";
        sout << " residual  " << std::setw(5) << std::setprecision(3) << resPull.first << " pull " << std::setw(6) << resPull.second << " "
             << m_idHelperSvc->toString(id);
        return sout.str();
    }

    template <class T> struct IdDataVec {
        typedef T Entry;
        typedef std::vector<Entry> EntryVec;
        typedef typename EntryVec::iterator EntryIt;

        IdDataVec() {}
        IdDataVec(const Identifier& i) : id(i) {}

        Identifier id;
        EntryVec data;
    };

    template <class T> struct SortIdDataVec {
        bool operator()(const IdDataVec<T>& d1, const IdDataVec<T>& d2) { return d1.id < d2.id; }
    };

    struct SortClByPull {
        bool operator()(const std::pair<double, DCMathSegmentMaker::Cluster2D>& d1,
                        const std::pair<double, DCMathSegmentMaker::Cluster2D>& d2) {
            return std::abs(d1.first) < std::abs(d2.first);
        }
    };

    std::pair<std::pair<int, int>, bool> DCMathSegmentMaker::associateClustersToSegment(
        const TrkDriftCircleMath::Segment& segment, const Identifier& chid, const Amg::Transform3D& gToStation, ClusterVecPair& spVecs,
        double phimin, double phimax, std::vector<std::pair<double, std::unique_ptr<const Trk::MeasurementBase>> >& rioDistVec) const {
        typedef IdDataVec<std::pair<double, Cluster2D> > GasGapData;
        typedef IdDataVec<GasGapData> ChamberData;
        typedef std::vector<ChamberData> ChamberDataVec;
        typedef ChamberDataVec::iterator ChamberDataIt;
        ChamberDataVec chamberDataVec;
        bool isEndcap = m_idHelperSvc->isEndcap(chid);

        // keep track of the number of eta/phi hits on the segment
        bool refit = false;
        std::pair<std::pair<int, int>, bool> netaPhiHits(std::pair<int, int>(0, 0), false);
        if (segment.clusters().empty()) return netaPhiHits;

        std::vector<const Trk::MeasurementBase*> phiHits;

        // only refit if there are sufficient phi hits and no multiple phi hits per surface
        refit = true;

        // keep track of detector elements which space points added to the track
        std::set<Identifier> detElOnSegments;
        std::set<MuonStationIndex::PhiIndex> phiIndices;

        ATH_MSG_DEBUG(" Associating clusters: " << segment.clusters().size() << " number of space points " << spVecs.first.size());

        // associate space points and sort them per detector element and gas gap
        for (TrkDriftCircleMath::CLCit ccit = segment.clusters().begin(); ccit != segment.clusters().end(); ++ccit) {
            ATH_MSG_VERBOSE(" accessing cluster: " << ccit->index());
            const Cluster2D& spacePoint = spVecs.first[ccit->index()];

            // skip corrupt space points
            if (spacePoint.corrupt()) {
                ATH_MSG_DEBUG(" Found corrupt space point: index " << ccit->index());
                continue;
            }
            // reject TGC clusters that are not 2D
            if (m_reject1DTgcSpacePoints && !spacePoint.is2D() && m_idHelperSvc->isTgc(spacePoint.identify())) {
                ATH_MSG_DEBUG(" Rejecting 1D tgc space point " << m_idHelperSvc->toString(spacePoint.identify()));
                continue;
            }
            if (m_assumePointingPhi && spacePoint.is2D() && !checkPhiConsistency(spacePoint.globalPos.phi(), phimin, phimax)) {
                ATH_MSG_DEBUG(" Inconsistent phi angle, dropping space point: phi " << spacePoint.globalPos.phi() << " range " << phimin
                                                                                    << " " << phimax);
                continue;
            }

            std::pair<double, double> resPull = residualAndPullWithSegment(segment, spacePoint, gToStation);

            // if empty or new chamber, add chamber
            if (chamberDataVec.empty() || chamberDataVec.back().id != spacePoint.detElId) {
                detElOnSegments.insert(spacePoint.detElId);
                chamberDataVec.push_back(ChamberData(spacePoint.detElId));
                MuonStationIndex::PhiIndex phiIndex = m_idHelperSvc->phiIndex(spacePoint.detElId);
                phiIndices.insert(phiIndex);
            }

            // reference to current chamber data
            ChamberData& chamber = chamberDataVec.back();

            // if same detector element
            if (spacePoint.detElId == chamber.id) {
                // if chamber empty or new gas gap, add gasp gap
                if (chamber.data.empty() || chamber.data.back().id != spacePoint.gasGapId) {
                    chamber.data.push_back(GasGapData(spacePoint.gasGapId));
                }
            }

            // reference to current gas gap data
            GasGapData& gasGap = chamber.data.back();
            gasGap.data.push_back(std::make_pair(resPull.second, spacePoint));
        }

        // calculate the distance between the first and last station, use r in barrel and z in endcaps
        double posFirstPhiStation = 1e6;
        double posLastPhiStation = 0;

        // loop over chambers and create competing ROTs per chamber
        ChamberDataIt cit = chamberDataVec.begin();
        ChamberDataIt cit_end = chamberDataVec.end();
        for (; cit != cit_end; ++cit) {
            // select best clusters per gas gap in chamber
            std::list<const Trk::PrepRawData*> etaClusterVec;
            std::list<const Trk::PrepRawData*> phiClusterVec;
            std::set<Identifier> etaIds;
            // loop over gas gaps
            ChamberData::EntryIt git = cit->data.begin();
            ChamberData::EntryIt git_end = cit->data.end();
            for (; git != git_end; ++git) {
                // sort space points by their pull with the segment
                std::sort(git->data.begin(), git->data.end(), SortClByPull());

                // select all space points with a pull that is within 1 of the best pull
                double bestPull = std::abs(git->data.front().first);

                // count number of associated clusters in gas gap
                unsigned int nassociatedSp = 0;
                GasGapData::EntryIt cl_it = git->data.begin();
                while (cl_it != git->data.end() && std::abs(cl_it->first) - bestPull < 1.) {
                    Cluster2D& sp = cl_it->second;

                    // calculate distance to segment
                    double dist = distanceToSegment(segment, sp.globalPos, gToStation);
                    ATH_MSG_VERBOSE("    selected space point:  " << m_idHelperSvc->toString(sp.identify()) << " pull "
                                                                  << std::abs(cl_it->first) << " distance to segment " << dist << " phi "
                                                                  << sp.globalPos.phi());

                    // here keep open the option not to create CompetingMuonClustersOnTrack
                    if (sp.etaHit) {
                        if (!etaIds.count(sp.etaHit->identify())) {
                            etaIds.insert(sp.etaHit->identify());

                            if (m_createCompetingROTsEta)
                                etaClusterVec.push_back(sp.etaHit->prepRawData());
                            else {
                                rioDistVec.push_back(std::make_pair(dist, sp.etaHit->uniqueClone()));
                                ++netaPhiHits.first.first;
                            }
                        }
                    }
                    if (!sp.phiHits.empty()) {
                        if (m_createCompetingROTsPhi) {
                            // can have multiple phi hits per cluster, loop over phi hits and add them
                            for (std::vector<const MuonClusterOnTrack*>::iterator pit = sp.phiHits.begin(); pit != sp.phiHits.end();
                                 ++pit) {
                                phiClusterVec.push_back((*pit)->prepRawData());
                            }
                        } else {
                            // can have multiple phi hits per cluster, loop over phi hits and add them
                            for (std::vector<const MuonClusterOnTrack*>::iterator pit = sp.phiHits.begin(); pit != sp.phiHits.end();
                                 ++pit) {
                                rioDistVec.push_back(std::make_pair(dist, (*pit)->uniqueClone()));
                                ++netaPhiHits.first.second;
                                phiHits.push_back(*pit);

                                // calculate position
                                double phiPos =
                                    isEndcap ? std::abs((*pit)->globalPosition().z()) : std::abs((*pit)->globalPosition().perp());
                                if (phiPos < posFirstPhiStation) posFirstPhiStation = phiPos;
                                if (phiPos > posLastPhiStation) posLastPhiStation = phiPos;
                            }
                            if (sp.phiHits.size() > 1) refit = false;
                        }
                    }
                    ++nassociatedSp;
                    ++cl_it;
                }
                // multiple clusters in same gas gap, don't refit
                if (!m_createCompetingROTsPhi && nassociatedSp > 1) refit = false;
            }

            if (m_createCompetingROTsEta) {
                // create competing ROT for eta hits
                if (!etaClusterVec.empty()) {
                    std::unique_ptr<const CompetingMuonClustersOnTrack> etaCompCluster = m_compClusterCreator->createBroadCluster(etaClusterVec, 0.);
                    if (!etaCompCluster) {
                        ATH_MSG_DEBUG(" failed to create competing ETA ROT " << etaClusterVec.size());
                    } else {
                        double dist = distanceToSegment(segment, etaCompCluster->globalPosition(), gToStation);                       
                        ++netaPhiHits.first.first;
                        if (msgLvl(MSG::VERBOSE)) {
                            ATH_MSG_VERBOSE("    selected cluster:  " << m_idHelperSvc->toString(etaClusterVec.front()->identify()));
                            for (unsigned int i = 0; i < etaCompCluster->containedROTs().size(); ++i) {
                                ATH_MSG_VERBOSE(
                                    "               content:  " << m_idHelperSvc->toString(etaCompCluster->containedROTs()[i]->identify()));
                            }
                        }
                        rioDistVec.push_back(std::make_pair(dist, std::move(etaCompCluster)));
                    }
                }
            }

            if (m_createCompetingROTsPhi) {
                // create competing ROT for phi hits
                if (!phiClusterVec.empty()) {
                    std::unique_ptr<const CompetingMuonClustersOnTrack> phiCompCluster = m_compClusterCreator->createBroadCluster(phiClusterVec, 0.);
                    if (!phiCompCluster) {
                        ATH_MSG_DEBUG(" failed to create competing PHI ROT " << phiClusterVec.size());
                    } else {
                        double dist = distanceToSegment(segment, phiCompCluster->globalPosition(), gToStation);
                        phiHits.push_back(phiCompCluster.get());
                        
                        ++netaPhiHits.first.second;

                        if (msgLvl(MSG::VERBOSE)) {
                            ATH_MSG_VERBOSE("    selected cluster:  " << m_idHelperSvc->toString(phiClusterVec.front()->identify()));
                            for (unsigned int i = 0; i < phiCompCluster->containedROTs().size(); ++i) {
                                ATH_MSG_VERBOSE(
                                    "               content:  " << m_idHelperSvc->toString(phiCompCluster->containedROTs()[i]->identify()));
                            }
                        }
                        

                        // calculate position
                        double phiPos =
                            isEndcap ? std::abs(phiCompCluster->globalPosition().z()) : phiCompCluster->globalPosition().perp();
                        posFirstPhiStation = std::min(phiPos,posFirstPhiStation);
                        posLastPhiStation = std::max(phiPos,posLastPhiStation);
                        rioDistVec.push_back(std::make_pair(dist, std::move(phiCompCluster)));
                      
                    }
                }
            }
        }

        // add phi hits that were not associated with an eta hit (only in barrel)
        if ((!spVecs.second.empty() || m_recoverBadRpcCabling) && m_addUnassociatedPhiHits && !isEndcap) {
            TrkDriftCircleMath::ResidualWithSegment resWithSegment(segment);

            std::map<Identifier, std::list<const Trk::PrepRawData*> > phiClusterMap;

            std::set<const MuonClusterOnTrack*> selectedClusters;
            std::vector<Cluster2D*> phiClusters;
            phiClusters.reserve(spVecs.second.size());

            // create lists of PrepRawData per detector element
            ClusterIt pit = spVecs.second.begin();
            ClusterIt pit_end = spVecs.second.end();
            for (; pit != pit_end; ++pit) {
                if (!pit->phiHit || pit->corrupt()) {
                    ATH_MSG_WARNING(" phi clusters without phi hit!!");
                    continue;
                }
                phiClusters.push_back(&*pit);
                selectedClusters.insert(pit->phiHit);
            }

            unsigned int recoveredUnassociatedPhiHits(0);
            if (m_recoverBadRpcCabling) {
                // now loop over 2D space points and add the phi hits to the list if the detEl is not yet added to the
                // segment
                ClusterIt spit = spVecs.first.begin();
                ClusterIt spit_end = spVecs.first.end();
                for (; spit != spit_end; ++spit) {
                    // skip clusters without phi hit
                    if (!spit->phiHit || spit->corrupt()) continue;

                    // skip clusters in detector elements that are already associated (ok as this is only done for RPC)
                    if (detElOnSegments.count(spit->detElId)) continue;

                    MuonStationIndex::PhiIndex phiIndex = m_idHelperSvc->phiIndex(spit->detElId);
                    // skip clusters in detector layer
                    if (phiIndices.count(phiIndex)) continue;

                    bool wasFound = false;
                    for (std::vector<const MuonClusterOnTrack*>::iterator it = spit->phiHits.begin(); it != spit->phiHits.end(); ++it) {
                        // now to avoid duplicate also skip if the given ROT is already in the list
                        std::pair<std::set<const MuonClusterOnTrack*>::iterator, bool> insertResult = selectedClusters.insert(*it);
                        if (!insertResult.second) {
                            // flag as found
                            wasFound = true;

                            // just because I'm paranoid, remove the hits from this cluster that were already inserted up to
                            // here
                            for (std::vector<const MuonClusterOnTrack*>::iterator rmit = spit->phiHits.begin(); rmit != it; ++rmit)
                                selectedClusters.erase(*rmit);
                            break;
                        }
                    }
                    if (wasFound) continue;

                    // if we get here we should add the hit
                    phiClusters.push_back(&*spit);
                    ++recoveredUnassociatedPhiHits;
                }
            }

            unsigned int addedPhiHits(0);
            std::vector<Cluster2D*>::iterator cit = phiClusters.begin();
            std::vector<Cluster2D*>::iterator cit_end = phiClusters.end();
            for (; cit != cit_end; ++cit) {
                const Identifier& detElId = (*cit)->detElId;

                // check that detector element is not already added to segment
                if (detElOnSegments.count(detElId)) continue;

                MuonStationIndex::PhiIndex phiIndex = m_idHelperSvc->phiIndex(detElId);
                // skip clusters in detector layer
                if (phiIndices.count(phiIndex)) continue;

                // calculate local cluster position
                Amg::Vector3D locPos = gToStation * (*cit)->globalPos;

                // calculate intersect of segment with cluster
                TrkDriftCircleMath::Cluster cl(TrkDriftCircleMath::LocVec2D(locPos.y(), locPos.z()), 1.);
                double residual = resWithSegment.residual(cl);
                double segError = sqrt(resWithSegment.trackError2(cl));
                const MuonGM::RpcReadoutElement* detEl = dynamic_cast<const MuonGM::RpcReadoutElement*>((*cit)->phiHit->detectorElement());
                if (!detEl) {
                    ATH_MSG_WARNING("found RPC prd without RpcReadoutElement");
                    continue;
                }

                // perform bound check
                double stripLength = detEl->StripLength(1);
                bool inBounds = std::abs(residual) < 0.5 * stripLength + 2. + segError ? true : false;
                if (msgLvl(MSG::DEBUG)) {
                    ATH_MSG_DEBUG(" Unassociated " << m_idHelperSvc->toString((*cit)->phiHit->identify()) << " pos x " << cl.position().x()
                                                   << " pos y " << cl.position().y() << " : residual " << residual << " strip half length "
                                                   << 0.5 * stripLength << " segment error " << segError);
                    if (inBounds)
                        ATH_MSG_DEBUG(" inBounds");
                    else
                        ATH_MSG_DEBUG(" outBounds");
                }
                if (inBounds) {
                    // can have multiple phi hits per cluster, loop over phi hits and add them
                    for (std::vector<const MuonClusterOnTrack*>::iterator pit = (*cit)->phiHits.begin(); pit != (*cit)->phiHits.end();
                         ++pit) {
                        phiClusterMap[detElId].push_back((*pit)->prepRawData());
                    }
                }
            }

            // loop over detector elements and created competing ROTs

            std::map<Identifier, std::list<const Trk::PrepRawData*> >::iterator chit = phiClusterMap.begin();
            std::map<Identifier, std::list<const Trk::PrepRawData*> >::iterator chit_end = phiClusterMap.end();
            for (; chit != chit_end; ++chit) {
                std::list<const Trk::PrepRawData*>& prds = chit->second;
                if (prds.empty()) {
                    ATH_MSG_WARNING(" chamber without phi hits " << m_idHelperSvc->toString(chit->first));
                    continue;
                }

                std::unique_ptr<const CompetingMuonClustersOnTrack> phiCompCluster = m_compClusterCreator->createBroadCluster(prds, 0.);
                if (!phiCompCluster) {
                    ATH_MSG_DEBUG(" failed to create competing PHI ROT " << prds.size());
                } else {
                    double dist = distanceToSegment(segment, phiCompCluster->globalPosition(), gToStation);

                    if (std::abs(dist) > m_maxAssociateClusterDistance) {
                        
                        ATH_MSG_VERBOSE("    rejected unassociated cluster:  " << m_idHelperSvc->toString(prds.front()->identify())
                                                                               << "  distance to segment " << dist);
                        continue;
                    }
                    phiHits.push_back(phiCompCluster.get());
                    ++netaPhiHits.first.second;
                    ++addedPhiHits;
                    if (msgLvl(MSG::VERBOSE)) {
                        ATH_MSG_VERBOSE("    selected unassociated cluster:  " << m_idHelperSvc->toString(prds.front()->identify())
                                                                               << "  distance to segment " << dist);
                        for (unsigned int i = 0; i < phiCompCluster->containedROTs().size(); ++i) {
                            ATH_MSG_VERBOSE(
                                "               content:  " << m_idHelperSvc->toString(phiCompCluster->containedROTs()[i]->identify()));
                        }
                    }
                    rioDistVec.push_back(std::make_pair(dist, std::move(phiCompCluster)));
                }
            }
            ATH_MSG_VERBOSE("Added " << addedPhiHits << " unass phi hits out of " << spVecs.second.size()
                                     << " phi hits without eta hit and " << recoveredUnassociatedPhiHits << " with unassociated eta hit ");
        }

        // calculate distance between first and last phi trigger hit, refit if there is a good constraint on phi
        double phiDistanceMax = posLastPhiStation - posFirstPhiStation;
        if (isEndcap && phiDistanceMax < 1000.)
            refit = false;
        else if (phiDistanceMax < 400.)
            refit = false;

        netaPhiHits.second = refit;
        return netaPhiHits;
    }

    double DCMathSegmentMaker::distanceToSegment(const TrkDriftCircleMath::Segment& segment, const Amg::Vector3D& hitPos,
                                                 const Amg::Transform3D& gToStation) const {
        const TrkDriftCircleMath::Line& line = segment.line();
        TrkDriftCircleMath::TransformToLine toLine(line);
        double cos_sinLine = cot(line.phi());

        // calculate local AMDB position
        Amg::Vector3D locPos = gToStation * hitPos;

        TrkDriftCircleMath::LocVec2D lpos(locPos.y(), locPos.z());

        // calculate distance of segment to measurement surface
        double delta_y = lpos.y() - line.position().y();

        // calculate position of hit in line frame
        TrkDriftCircleMath::LocVec2D lineSurfaceIntersect(delta_y * cos_sinLine + line.position().x(), lpos.y());

        // calculate position of hit in line frame
        TrkDriftCircleMath::LocVec2D pointOnHit = toLine.toLine(lineSurfaceIntersect);

        return pointOnHit.x();
    }

    std::unique_ptr<DataVector<const Trk::MeasurementBase>> DCMathSegmentMaker::createROTVec(
        std::vector<std::pair<double,  std::unique_ptr<const Trk::MeasurementBase>> >& rioDistVec) const {
        // sort hits according to they distance to the segment position
        std::sort(rioDistVec.begin(), rioDistVec.end(), SortByDistanceToSegment());

        std::unique_ptr<DataVector<const Trk::MeasurementBase>> rioVec = std::make_unique<DataVector<const Trk::MeasurementBase>>(SG::OWN_ELEMENTS);
        rioVec->reserve(rioDistVec.size());
        for (std::pair<double, std::unique_ptr<const Trk::MeasurementBase>>& rdit : rioDistVec) { rioVec->push_back(std::move(rdit.second)); }
        rioDistVec.clear();
        return rioVec;
    }

    std::pair<double, double> DCMathSegmentMaker::residualAndPullWithSegment(const TrkDriftCircleMath::Segment& segment,
                                                                             const Cluster2D& spacePoint,
                                                                             const Amg::Transform3D& gToStation) const {
        const TrkDriftCircleMath::Line& line = segment.line();
        double cos_sinLine = cot(line.phi());

        // calculate sp postion in AMDB reference frame
        Amg::Vector3D locPos = gToStation * spacePoint.globalPos;
        TrkDriftCircleMath::LocVec2D lpos(locPos.y(), locPos.z());

        // calculate distance of segment to measurement surface
        double delta_y = lpos.y() - line.position().y();

        // calculate position of hit in line frame
        TrkDriftCircleMath::LocVec2D lineSurfaceIntersect(delta_y * cos_sinLine + line.position().x(), lpos.y());

        // calculate position of hit in line frame
        double residual = lpos.x() - lineSurfaceIntersect.x();
        double pull = residual / spacePoint.error;
        return std::make_pair(residual, pull);
    }

    std::vector<Identifier> DCMathSegmentMaker::calculateHoles(const EventContext& ctx,
        Identifier chid, const Amg::Vector3D& gpos, const Amg::Vector3D& gdir, bool hasMeasuredCoordinate, std::set<Identifier>& deltaVec,
        std::set<Identifier>& outoftimeVec, const  std::vector<std::pair<double,  std::unique_ptr<const Trk::MeasurementBase>> >& rioDistVec) const {
        // calculate crossed tubes
        SG::ReadCondHandle<Muon::MuonIntersectGeoData> InterSectSvc{m_chamberGeoKey, ctx};
        if (!InterSectSvc.isValid()) {
            ATH_MSG_ERROR("Null pointer to the read MuonDetectorManager conditions object");
            return {};
        }
        const MuonStationIntersect intersect = InterSectSvc->tubesCrossedByTrack(chid, gpos, gdir);
        const MuonGM::MuonDetectorManager* MuonDetMgr = InterSectSvc->detMgr();

        // set to identify the hit on the segment
        std::set<Identifier> hitsOnSegment, chambersOnSegment;
        int firstLayer{-1}, lastLayer{-1};
        for (const std::pair<double, std::unique_ptr<const Trk::MeasurementBase>>& rdit : rioDistVec) {
            const MdtDriftCircleOnTrack* mdt = dynamic_cast<const MdtDriftCircleOnTrack*>(rdit.second.get());
            if (mdt) {
                const Identifier& id = mdt->identify();
                int layer = (m_idHelperSvc->mdtIdHelper().tubeLayer(id) - 1) + 4 * (m_idHelperSvc->mdtIdHelper().multilayer(id) - 1);
                if (firstLayer == -1)
                    firstLayer = layer;
                else
                    lastLayer = layer;

                hitsOnSegment.insert(mdt->identify());
                chambersOnSegment.insert(m_idHelperSvc->chamberId(mdt->identify()));
            }
        }

        // cross check for cosmic case
        if (firstLayer > lastLayer) { std::swap(firstLayer, lastLayer); }
        ATH_MSG_VERBOSE(" Tube layer ranges: " << firstLayer << " -- " << lastLayer << " crossed tubes "
                                               << intersect.tubeIntersects().size());
        // clear hole vector
        std::vector<Identifier> holeVec;
        for (unsigned int ii = 0; ii < intersect.tubeIntersects().size(); ++ii) {
            const MuonTubeIntersect& tint = intersect.tubeIntersects()[ii];
            if (!chambersOnSegment.count(m_idHelperSvc->chamberId(tint.tubeId))) {
                ATH_MSG_VERBOSE(" chamber not on segment, not counting tube  " << tint.rIntersect << " l " << tint.xIntersect << " "
                                                                               << m_idHelperSvc->toString(tint.tubeId));
                continue;
            }

            const Identifier& id = tint.tubeId;
            int layer = (m_idHelperSvc->mdtIdHelper().tubeLayer(id) - 1) + 4 * (m_idHelperSvc->mdtIdHelper().multilayer(id) - 1);

            bool notBetweenHits = layer < firstLayer || layer > lastLayer;
            double distanceCut = hasMeasuredCoordinate ? -20 : -200.;
            double innerRadius = MuonDetMgr->getMdtReadoutElement(id)->innerTubeRadius();
            if (notBetweenHits && (std::abs(tint.rIntersect) > innerRadius || (!m_allMdtHoles && tint.xIntersect > distanceCut))) {
                ATH_MSG_VERBOSE(" not counting tube:  distance to wire " << tint.rIntersect << " dist to tube end " << tint.xIntersect
                                                                         << " " << m_idHelperSvc->toString(tint.tubeId));
            } else {
                // check whether there is a hit in this tube
                if (hitsOnSegment.count(tint.tubeId)) {
                    ATH_MSG_VERBOSE(" tube on segment:  distance to wire " << tint.rIntersect << " dist to tube end " << tint.xIntersect
                                                                           << " " << m_idHelperSvc->toString(tint.tubeId));
                    continue;
                }
                // check whether there is a delta electron in this tube
                if (m_removeDeltas) {
                    if (deltaVec.count(tint.tubeId)) {
                        ATH_MSG_VERBOSE(" removing delta, distance to wire " << tint.rIntersect << " dist to tube end " << tint.xIntersect
                                                                             << " " << m_idHelperSvc->toString(tint.tubeId));
                        continue;
                    }

                    const MdtPrepData* prd = findMdt(id);
                    if (prd && std::abs(prd->localPosition()[Trk::locR]) < std::abs(tint.rIntersect)) {
                        ATH_MSG_VERBOSE(" found and removed delta, distance to wire " << tint.rIntersect << " dist to tube end "
                                                                                      << tint.xIntersect << " "
                                                                                      << m_idHelperSvc->toString(tint.tubeId));
                        continue;
                    }
                }
                if (msgLvl(MSG::VERBOSE)) {
                    std::string holeType = "hole ";
                    if (outoftimeVec.count(tint.tubeId)) holeType = "Out-of-time ";
                    ATH_MSG_VERBOSE(holeType << " distance to wire " << tint.rIntersect << " dist to tube end " << tint.xIntersect << " "
                                             << m_idHelperSvc->toString(tint.tubeId));
                    if (!notBetweenHits) ATH_MSG_VERBOSE(" between hits");
                }

                holeVec.push_back(tint.tubeId);
            }
        }
        return holeVec;
    }

    const MdtPrepData* DCMathSegmentMaker::findMdt(const Identifier& id) const {
        IdentifierHash colHash;
        if (m_idHelperSvc->mdtIdHelper().get_module_hash(m_idHelperSvc->chamberId(id), colHash) == 0) {
            const MdtPrepDataContainer* prdContainer = nullptr;
            SG::ReadHandle<Muon::MdtPrepDataContainer> MdtCont(m_mdtKey);
            if (!MdtCont.isValid()) {
                ATH_MSG_WARNING("Cannot retrieve MdtPrepDataContainer ");
                return nullptr;
            } else {
                prdContainer = MdtCont.cptr();
            }
            if (!prdContainer) {
                ATH_MSG_WARNING("No MDT prd collection retrieved");
                return nullptr;
            }
            auto collptr = prdContainer->indexFindPtr(colHash);
            if (collptr == nullptr) return nullptr;
            const MdtPrepDataCollection& col = *collptr;
            MdtPrepDataCollection::const_iterator mdtIt = col.begin();
            MdtPrepDataCollection::const_iterator mdtIt_end = col.end();
            for (; mdtIt != mdtIt_end; ++mdtIt) {
                if ((*mdtIt)->identify() == id) return *mdtIt;
            }
        }
        return nullptr;
    }

    const MdtDriftCircleOnTrack* DCMathSegmentMaker::findFirstRotInChamberWithMostHits(
        const std::vector<const MdtDriftCircleOnTrack*>& mdts) const {
        int hitsInChamberWithMostHits = 0;
        std::map<Identifier, int> hitsPerChamber;
        int currentSector = -1;
        const MdtDriftCircleOnTrack* rotInChamberWithMostHits = nullptr;

        // loop over all MDTs and count number of MDTs per chamber
        for (std::vector<const MdtDriftCircleOnTrack*>::const_iterator it = mdts.begin(); it != mdts.end(); ++it) {
            const MdtDriftCircleOnTrack* rot = *it;

            if (!rot) {
                ATH_MSG_WARNING(" rot not a MdtDriftCircleOnTrack ");
                continue;
            }
            Identifier chId = m_idHelperSvc->chamberId((*it)->identify());
            int sector = m_idHelperSvc->sector(chId);
            if (currentSector == -1) {
                currentSector = sector;
            } else if (sector != currentSector) {
                return nullptr;
            }
            int& hitsInCh = hitsPerChamber[chId];
            ++hitsInCh;
            if (hitsInCh > hitsInChamberWithMostHits) {
                hitsInChamberWithMostHits = hitsInCh;
                rotInChamberWithMostHits = rot;
            }
        }
        return rotInChamberWithMostHits;
    }

    std::pair<bool, double> DCMathSegmentMaker::rotateLocalAngleXZIntoBounds(double xline, double zline, double dXdZ,
                                                                             std::vector<DCMathSegmentMaker::HitInXZ>& hits) const {
        // routine checks whether hits are in bounds, if not rotates the local angle XZ until all hits are in bounds

        std::pair<bool, double> result(false, dXdZ);
        if (hits.empty()) return result;

        double dXdZMinLeft{1e9}, dXdZMinRight{1e9};
        std::vector<HitInXZ>::iterator hit = hits.begin();
        std::vector<HitInXZ>::iterator hit_end = hits.end();
        ATH_MSG_VERBOSE(" recalculating using phi measurement " << m_idHelperSvc->toString(hits.front().id) << " pos (" << xline << ","
                                                                << zline << ")");
        for (; hit != hit_end; ++hit) {
            // sl extrapolation
            double locExX = xline + dXdZ * (hit->z - zline);

            bool newBoundLeft = false;
            bool newBoundRight = false;

            // sanity checks to avoid division + only for MDT for now
            if (hit->isMdt) {
                // now calculate, if out of bounds, for which angle the situation would be ok
                double dz = (hit->z - zline);
                if (std::abs(dz) > 0.01) {
                    if (locExX < hit->xmin) {
                        double inBoundsdXdZ = (hit->xmin - xline) / dz;
                        if (std::abs(inBoundsdXdZ) < std::abs(dXdZMinLeft)) dXdZMinLeft = inBoundsdXdZ;
                        newBoundLeft = true;
                    }
                    if (locExX > hit->xmax) {
                        double inBoundsdXdZ = (hit->xmax - xline) / dz;
                        if (std::abs(inBoundsdXdZ) < std::abs(dXdZMinRight)) dXdZMinRight = inBoundsdXdZ;
                        newBoundRight = true;
                    }
                }
            }
            if (msgLvl(MSG::VERBOSE)) {
                ATH_MSG_VERBOSE("  " << std::setw(65) << m_idHelperSvc->toString(hit->id) << " pos (" << std::setw(6) << (int)hit->x << ","
                                     << std::setw(6) << (int)hit->z << ")  ex pos " << std::setw(6) << (int)locExX << " min "
                                     << std::setw(6) << (int)hit->xmin << " max " << std::setw(6) << (int)hit->xmax << " phimin "
                                     << std::setw(6) << hit->phimin << " phimax " << std::setw(6) << hit->phimax);
                if (locExX < hit->xmin || locExX > hit->xmax) ATH_MSG_VERBOSE(" outBounds ");
                if (newBoundRight) ATH_MSG_VERBOSE(" right " << dXdZMinRight);
                if (newBoundLeft) ATH_MSG_VERBOSE(" left " << dXdZMinLeft);
            }
        }

        if (msgLvl(MSG::VERBOSE)) {
            if (std::abs(dXdZMinRight) < 1e8) ATH_MSG_VERBOSE(" selected right bound " << dXdZMinRight);
            if (std::abs(dXdZMinLeft) < 1e8) ATH_MSG_VERBOSE(" selected left bound " << dXdZMinLeft);
        }

        if (std::abs(dXdZMinRight) < 1e8) {
            result.first = true;
            result.second = dXdZMinRight;
        } else if (std::abs(dXdZMinLeft) < 1e8) {
            result.first = true;
            result.second = dXdZMinLeft;
        }

        return result;
    }

    bool DCMathSegmentMaker::checkBoundsInXZ(double xline, double zline, double dXdZ,
                                             std::vector<DCMathSegmentMaker::HitInXZ>& hits) const {
        bool ok = true;

        // look over hits and check whether all are in bounds
        std::vector<HitInXZ>::iterator hit = hits.begin();
        std::vector<HitInXZ>::iterator hit_end = hits.end();
        for (; hit != hit_end; ++hit) {
            bool outBounds = false;
            double locExX = xline + dXdZ * (hit->z - zline);
            if (hit->isMdt && (locExX < hit->xmin - 1. || locExX > hit->xmax + 1.)) {
                ok = false;
                outBounds = true;
                if (!msgLvl(MSG::DEBUG)) break;
            }

            if (outBounds && msgLvl(MSG::DEBUG)) {
                ATH_MSG_DEBUG("  " << std::setw(65) << m_idHelperSvc->toString(hit->id) << " pos (" << std::setw(6) << (int)hit->x << ","
                                   << std::setw(6) << (int)hit->z << ")  ex pos " << std::setw(6) << (int)locExX << " min " << std::setw(6)
                                   << (int)hit->xmin << " max " << std::setw(6) << (int)hit->xmax << " phimin " << std::setw(6)
                                   << hit->phimin << " phimax " << std::setw(6) << hit->phimax << " outBounds, cross-check");
            }
        }
        return ok;
    }

    bool DCMathSegmentMaker::updateSegmentPhi(const Amg::Vector3D& gpos, const Amg::Vector3D& gdir, Amg::Vector2D& segLocPos,
                                              Trk::LocalDirection& segLocDir, Trk::PlaneSurface& surf,
                                              const std::vector<const Trk::MeasurementBase*>& rots, double seg_phimin,
                                              double seg_phimax) const {
        bool hasUpdated = false;

        const Amg::Transform3D& segmentToGlobal = surf.transform();
        Amg::Transform3D gToSegment = surf.transform().inverse();
        Amg::Vector3D ldir = gToSegment * gdir;

        // ensure that we can extrapolate
        if (ldir.z() < 0.0001) return false;

        double dXdZ = ldir.x() / ldir.z();
        double dYdZ = ldir.y() / ldir.z();
        Amg::Vector3D lsegPos = gToSegment * gpos;
        double xline = lsegPos.x();
        double yline = lsegPos.y();
        double zline = lsegPos.z();
        ATH_MSG_VERBOSE(" Associated  hits " << rots.size() << " angleXZ " << 90. * segLocDir.angleXZ() / acos(0.) << " dXdZ " << dXdZ
                                             << " seg Pos (" << xline << " " << zline << ") " << segLocPos);

        std::vector<HitInXZ> hits;
        hits.reserve(rots.size());

        unsigned int nphiHits(0);
        HitInXZ* firstPhiHit = nullptr;
        HitInXZ* lastPhiHit = nullptr;

        for (std::vector<const Trk::MeasurementBase*>::const_iterator it = rots.begin(); it != rots.end(); ++it) {
            Identifier id = m_edmHelperSvc->getIdentifier(**it);
            if (!id.is_valid()) continue;
            Amg::Vector3D lpos;
            double lxmin{0}, lxmax{0}, phimin{0}, phimax{0};
            bool isMdt = m_idHelperSvc->isMdt(id);
            bool measuresPhi = m_idHelperSvc->measuresPhi(id);
            if (isMdt) {
                lpos.setZero();
                const MdtDriftCircleOnTrack* mdt = dynamic_cast<const MdtDriftCircleOnTrack*>(*it);
                if (!mdt) continue;
                TubeEnds tubeEnds = localTubeEnds(*mdt, gToSegment, segmentToGlobal);

                lxmin = tubeEnds.lxmin;
                lxmax = tubeEnds.lxmax;
                phimin = tubeEnds.phimin;
                phimax = tubeEnds.phimax;
            } else {
                lpos = gToSegment * (*it)->globalPosition();
                lxmin = lpos.x() - 5 * Amg::error((*it)->localCovariance(), Trk::locX);
                lxmax = lpos.x() + 5 * Amg::error((*it)->localCovariance(), Trk::locX);

                const CompetingMuonClustersOnTrack* crot = dynamic_cast<const CompetingMuonClustersOnTrack*>(*it);
                if (!measuresPhi) {
                    if (crot) {
                        const MuonGM::RpcReadoutElement* detEl =
                            dynamic_cast<const MuonGM::RpcReadoutElement*>(crot->containedROTs().front()->prepRawData()->detectorElement());
                        if (detEl) {
                            // perform bound check
                            double stripLength = detEl->StripLength(0);

                            lxmin = lpos.x() - 0.5 * stripLength;
                            lxmax = lpos.x() + 0.5 * stripLength;
                        }
                    }
                    Amg::Vector3D locPosition = lpos;
                    locPosition[0] = lxmin;
                    Amg::Vector3D globalPos = segmentToGlobal * locPosition;
                    double phi1 = globalPos.phi();

                    locPosition[0] = lxmax;
                    globalPos = segmentToGlobal * locPosition;
                    double phi2 = globalPos.phi();
                    phimin = phi1 < phi2 ? phi1 : phi2;
                    phimax = phi1 < phi2 ? phi2 : phi1;

                } else {
                    if (m_idHelperSvc->isTgc(id)) {
                        // need some special tricks for TGC phi hits as their reference plane can be rotated
                        // with respect to the MDT frame

                        // get orientation angle of strip to rotate back from local frame to strip
                        // copy code from ROT creator
                        int stripNo = m_idHelperSvc->tgcIdHelper().channel(id);
                        int gasGap = m_idHelperSvc->tgcIdHelper().gasGap(id);
                        if (!crot) {
                            ATH_MSG_WARNING("dynamic cast failed for CompetingMuonClustersOnTrack");
                            continue;
                        }
                        const MuonGM::TgcReadoutElement* detEl =
                            dynamic_cast<const MuonGM::TgcReadoutElement*>(crot->containedROTs().front()->prepRawData()->detectorElement());
                        if (!detEl) {
                            ATH_MSG_WARNING("dynamic cast failed for TgcReadoutElement");
                            continue;
                        }
                        // calculate two points along the tgc phi strip in the local tgc reference frame
                        Amg::Vector3D lposTGC = detEl->localChannelPos(id);
                        double z_shift = lposTGC.z() + 10;
                        double locy_shift = detEl->stripCtrX(gasGap, stripNo, z_shift);
                        if (0 < detEl->getStationEta()) { locy_shift *= -1.; }
                        Amg::Vector3D lpos_shift(lposTGC.x(), locy_shift, z_shift);

                        // transform the two points to global coordinates
                        const Amg::Transform3D tgcTrans = detEl->absTransform();
                        Amg::Vector3D gposL = tgcTrans * lposTGC;
                        Amg::Vector3D gposL_shift = tgcTrans * lpos_shift;

                        // now transform them into the segment frame
                        Amg::Vector3D lposSeg = gToSegment * gposL;
                        Amg::Vector3D lposSeg_shift = gToSegment * gposL_shift;

                        // calculate the y coordinate of the intersect of the segment with the TGC plane in the segment
                        // frame
                        double segYAtHit = yline + dYdZ * (lposSeg.z() - zline);

                        // the TGC phi strip is a line in the xy plane, calculate the x position of the point on the line
                        // at the y intersect position of the segment
                        double tgcdX = lposSeg_shift.x() - lposSeg.x();
                        double tgcdY = lposSeg_shift.y() - lposSeg.y();
                        if (std::abs(tgcdY) < 0.0001) {
                            ATH_MSG_WARNING(" Bad TGC phi strip orientation ");
                            continue;
                        }
                        double tgcExX = tgcdX / tgcdY * (segYAtHit - lposSeg.y()) + lposSeg.x();
                        lpos[0] = tgcExX;
                        lpos[1] = segYAtHit;
                        if (msgLvl(MSG::VERBOSE))
                            ATH_MSG_VERBOSE(" In seg frame: phi pos " << lposSeg << " shifted pos " << lposSeg_shift
                                                                      << " intersect with segment " << lpos);
                    }
                    Amg::Vector3D globalPos = segmentToGlobal * lpos;
                    phimin = globalPos.phi();
                    phimax = phimin;

                    // check whether phi is consistent with segment phi range
                    bool phiOk = checkPhiConsistency(phimin, seg_phimin, seg_phimax);
                    if (!phiOk) {
                        if (msgLvl(MSG::DEBUG))
                            ATH_MSG_DEBUG(" Inconsistent phi " << phimin << " range " << seg_phimin << " " << seg_phimax);
                    }
                }
            }

            hits.push_back(HitInXZ(id, isMdt, measuresPhi, lpos.x(), lpos.z(), lxmin, lxmax, phimin, phimax));
            if (measuresPhi) {
                ++nphiHits;
                if (!firstPhiHit)
                    firstPhiHit = new HitInXZ(hits.back());
                else {
                    double distPhiHits = std::abs(firstPhiHit->z - hits.back().z);
                    if (distPhiHits > 500.) {
                        if (lastPhiHit) delete lastPhiHit;
                        lastPhiHit = new HitInXZ(hits.back());
                    } else {
                        // not count this phi hit
                        --nphiHits;
                        if (msgLvl(MSG::DEBUG)) { ATH_MSG_DEBUG(" close phi hits, distance " << distPhiHits); }
                    }
                }
            }
            if (msgLvl(MSG::VERBOSE)) {
                double locExX = xline + dXdZ * (lpos.z() - zline);
                ATH_MSG_VERBOSE("  " << std::setw(65) << m_idHelperSvc->toString(id) << " pos (" << std::setw(6) << (int)lpos.x() << ","
                                     << std::setw(6) << (int)lpos.z() << ")  ex pos " << std::setw(6) << (int)locExX << " min "
                                     << std::setw(6) << (int)lxmin << " max " << std::setw(6) << (int)lxmax << " phimin " << std::setw(6)
                                     << phimin << " phimax " << std::setw(6) << phimax);
                if (lpos.x() < lxmin || lpos.x() > lxmax) ATH_MSG_VERBOSE(" outBounds");
            }
        }

        if (nphiHits == 1) {
            if (!firstPhiHit) {
                ATH_MSG_WARNING(" Pointer to first phi hit not set, this should not happen! ");
            } else {
                if (xline != firstPhiHit->x) {
                    hasUpdated = true;

                    // use phi position of the phi hit
                    xline = firstPhiHit->x;
                    zline = firstPhiHit->z;

                    if (m_assumePointingPhi) {
                        Amg::Vector3D ipLocPos = gToSegment * Amg::Vector3D(0., 0., 0.);
                        if (msgLvl(MSG::VERBOSE)) ATH_MSG_VERBOSE(" IP position in local frame " << ipLocPos);

                        double dz = ipLocPos.z() - zline;
                        if (std::abs(dz) > 0.001) {
                            if (msgLvl(MSG::VERBOSE))
                                ATH_MSG_VERBOSE(" hit (" << xline << "," << zline << ")  IP (" << ipLocPos.x() << "," << ipLocPos.z()
                                                         << ")  dXdZ " << (ipLocPos.x() - xline) / dz << " old " << dXdZ);
                            dXdZ = (ipLocPos.x() - xline) / dz;
                        }
                    }
                }
            }
        } else if (nphiHits == 2) {
            if (!firstPhiHit || !lastPhiHit) {
                ATH_MSG_WARNING(" Pointer to one of the two phi hit not set, this should not happen! ");
            } else {
                double dz = lastPhiHit->z - firstPhiHit->z;
                // use phi position of the first hit
                xline = firstPhiHit->x;
                zline = firstPhiHit->z;
                if (std::abs(dz) > 300.) {
                    double dx = lastPhiHit->x - firstPhiHit->x;
                    hasUpdated = true;

                    // if the two hits are far enough apart, also use the direction of the line connecting the two hits.
                    dXdZ = dx / dz;
                }
            }
        } else {
            // in all other cases just rotate until the MDTs are ok
        }

        if (hasUpdated) {
            // move segment to position of phi hit
            double segX = xline - dXdZ * zline;

            // finally check whether now everything is in bounds
            bool ok = checkBoundsInXZ(segX, 0., dXdZ, hits);
            if (!ok) {
                // give WARNING and give up for now
                ATH_MSG_DEBUG("still several out of bounds hits after rotation: posx(" << segX << ") dXdZ " << dXdZ
                                                                                       << " keeping old result ");
            }

            // update segment parameters
            double alphaYZ = segLocDir.angleYZ();
            double alphaXZ = atan2(1, dXdZ);

            segLocPos[Trk::locX] = segX;
            segLocDir = Trk::LocalDirection(alphaXZ, alphaYZ);
        }

        delete firstPhiHit;
        delete lastPhiHit;
        return hasUpdated;
    }

    DCMathSegmentMaker::TubeEnds DCMathSegmentMaker::localTubeEnds(const MdtDriftCircleOnTrack& mdt, const Amg::Transform3D& gToSegment,
                                                                   const Amg::Transform3D& segmentToG) const {
        TubeEnds tubeEnds;
        const Identifier& id = mdt.identify();
        Amg::Vector3D lpos = gToSegment * mdt.prepRawData()->globalPosition();

        // use readout and hv side as the surface frame is not that of the chamber
        Amg::Vector3D lropos = gToSegment * mdt.prepRawData()->detectorElement()->ROPos(id);
        Amg::Vector3D lhvpos = lpos + (lpos - lropos);

        // rescale to correctly take into account active tube length
        double tubeLen = (lropos - lhvpos).mag();
        double activeTubeLen =
            mdt.detectorElement()->getActiveTubeLength(m_idHelperSvc->mdtIdHelper().tubeLayer(id), m_idHelperSvc->mdtIdHelper().tube(id));
        double scaleFactor = activeTubeLen / tubeLen;
        lropos[0] = scaleFactor * lropos.x();
        lhvpos[0] = scaleFactor * lhvpos.x();

        if (lropos.x() < lhvpos.x()) {
            tubeEnds.lxmin = lropos.x();
            tubeEnds.lxmax = lhvpos.x();
        } else {
            tubeEnds.lxmin = lhvpos.x();
            tubeEnds.lxmax = lropos.x();
        }

        Amg::Vector3D ropos = segmentToG * lropos;
        Amg::Vector3D hvpos = segmentToG * lhvpos;
        double phiRO = ropos.phi();
        double phiHV = hvpos.phi();
        tubeEnds.phimin = phiRO < phiHV ? phiRO : phiHV;
        tubeEnds.phimax = phiRO < phiHV ? phiHV : phiRO;
        return tubeEnds;
    }

    void DCMathSegmentMaker::updatePhiRanges(double phiminNew, double phimaxNew, double& phiminRef, double& phimaxRef) const {
        // check whether we are at the boundary where phi changes sign
        if (phiminRef * phimaxRef < 0.) {
            if (phiminRef < -1.1) {
                if (phiminRef > phiminNew) phiminRef = phiminNew;
                if (phimaxRef < phimaxNew) phimaxRef = phimaxNew;
            } else {
                if (phiminRef < phiminNew) phiminRef = phiminNew;
                if (phimaxRef > phimaxNew) phimaxRef = phimaxNew;
            }
        } else {
            // if not life is 'easy'
            if (phiminRef < 0.) {
                if (phiminRef < phiminNew) phiminRef = phiminNew;
                if (phimaxRef > phimaxNew) phimaxRef = phimaxNew;
            } else {
                if (phiminRef > phiminNew) phiminRef = phiminNew;
                if (phimaxRef < phimaxNew) phimaxRef = phimaxNew;
            }
        }
    }

    bool DCMathSegmentMaker::checkPhiConsistency(double phi, double phimin, double phimax) const {
        // only if assuming pointing phi
        if (!m_assumePointingPhi) return true;

        bool phiOk = true;
        double offset = 0.05;
        if (phimin * phimax < 0.) {
            if (phi < 0.) {
                if (phi > -1.1) {
                    if (phi < phimin - offset) phiOk = false;
                } else {
                    if (phi > phimin + offset) phiOk = false;
                }
            } else {
                if (phi > 1.1) {
                    if (phi < phimax - offset) phiOk = false;
                } else {
                    if (phi > phimax + offset) phiOk = false;
                }
            }
        } else {
            if (phi < phimin - offset || phi > phimax + offset) phiOk = false;
        }
        return phiOk;
    }

    Amg::Vector3D DCMathSegmentMaker::updateDirection(double linephi, const Trk::PlaneSurface& surf, const Amg::Vector3D& roaddir,
                                                      bool isCurvedSegment) const {
        // Local direction along precision measurement (0, dy, dz)
        Trk::LocalDirection segLocDirs(M_PI_2, linephi);
        Amg::Vector3D gdirs;
        surf.localToGlobalDirection(segLocDirs, gdirs);
        // Local direction in plane  (1,0,0)
        Trk::LocalDirection segLocDiro(0., M_PI_2);
        Amg::Vector3D gdiro;
        surf.localToGlobalDirection(segLocDiro, gdiro);

        // recalculate the value of the local XZ angle for the give YZ angle of the segment such that the global phi
        // direction remains unchanged
        double dx = std::sin(gdirs.theta()) * std::cos(gdirs.phi());
        double dy = std::sin(gdirs.theta()) * std::sin(gdirs.phi());
        double dz = std::cos(gdirs.theta());

        // vector gdiro

        double dxo = std::sin(gdiro.theta()) * std::cos(gdiro.phi());
        double dyo = std::sin(gdiro.theta()) * std::sin(gdiro.phi());
        double dzo = std::cos(gdiro.theta());

        // solve system real direction = A * gdir + B * gdiro
        //            phi global  constraint: (1)*sin(phi road) - (2)*cos(phi road) = 0
        //            ( A * dx + B * dxo ) sin (phi ) - (A * dy + B *dyo ) cos (phi) = 0
        //              A ( dx sin - dy cos ) + B (dx0 sin -dy0 cos) = A * a0 + B * b0 = 0
        //              psi = atan (-b0 , a0)

        double a0 = dx * std::sin(roaddir.phi()) - dy * std::cos(roaddir.phi());
        double b0 = dxo * std::sin(roaddir.phi()) - dyo * std::cos(roaddir.phi());
        if (b0 < 1e-8 && b0 > 0) b0 = 1e-8;
        if (b0 > -1e-8 && b0 < 0) b0 = -1e-8;
        double dxn = dx - a0 * dxo / b0;
        double dyn = dy - a0 * dyo / b0;
        double dzn = dz - a0 * dzo / b0;
        double norm = std::sqrt(dxn * dxn + dyn * dyn + dzn * dzn);

        // flip the sign if the direction NOT parallel to road
        if (m_assumePointingPhi) {
            if (dxn * roaddir.x() + dyn * roaddir.y() + dzn * roaddir.z() < 0.) { norm = -norm; }
        } else {
            if (dxn * roaddir.x() + dyn * roaddir.y() < 0.) { norm = -norm; }
        }

        if (isCurvedSegment) norm = norm / 2.;

        //
        // Follow segment fit direction
        //
        dxn = dxn / norm;
        dyn = dyn / norm;
        dzn = dzn / norm;

        return Amg::Vector3D(dxn, dyn, dzn);
    }
}  // namespace Muon
