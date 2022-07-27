/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonTrackToSegmentTool.h"

#include <set>

#include "EventPrimitives/EventPrimitivesHelpers.h"
#include "MuonRIO_OnTrack/CscClusterOnTrack.h"
#include "MuonRIO_OnTrack/MMClusterOnTrack.h"
#include "MuonRIO_OnTrack/MdtDriftCircleOnTrack.h"
#include "MuonReadoutGeometry/MdtReadoutElement.h"
#include "MuonSegment/MuonSegment.h"
#include "MuonSegment/MuonSegmentQuality.h"
#include "TrkEventPrimitives/JacobianPhiThetaLocalAngles.h"
#include "TrkEventPrimitives/LocalDirection.h"
#include "TrkGeometry/MagneticFieldProperties.h"
#include "TrkTrack/Track.h"

namespace Muon {
    MuonTrackToSegmentTool::MuonTrackToSegmentTool(const std::string& t, const std::string& n, const IInterface* p) : AthAlgTool(t, n, p) {
        declareInterface<IMuonTrackToSegmentTool>(this);
    }

    StatusCode MuonTrackToSegmentTool::initialize() {
        ATH_CHECK(m_propagator.retrieve());
        ATH_CHECK(m_idHelperSvc.retrieve());
        ATH_CHECK(m_edmHelperSvc.retrieve());
        ATH_CHECK(m_printer.retrieve());
        ATH_CHECK(m_chamberGeoKey.initialize());
        return StatusCode::SUCCESS;
    }

    MuonSegment* MuonTrackToSegmentTool::convert(const EventContext& ctx, const Trk::Track& track) const {
        /** convert track to segment, express the new segment parameters on the surface of the first segment */

        ATH_MSG_DEBUG(" creating MuonSegment from track ");

        const Trk::Perigee* perigee = track.perigeeParameters();
        if (!perigee) {
            ATH_MSG_WARNING(" was expecting a perigee here... ");
            return nullptr;
        }

        const Trk::FitQuality* fq = track.fitQuality();
        if (!fq) {
            ATH_MSG_WARNING(" was expecting a FitQuality here... ");
            return nullptr;
        }

        std::set<Identifier> chIds;

        // copy rots, get surface
        std::unique_ptr<DataVector<const Trk::MeasurementBase>> rots = std::make_unique<DataVector<const Trk::MeasurementBase>>();
        rots->reserve(track.measurementsOnTrack()->size());

        // loop over TSOS
        const Trk::TrackStates* states = track.trackStateOnSurfaces();
        if (!states) {
            ATH_MSG_WARNING(" track without states, discarding track ");
            return nullptr;
        }
        // track direction vector
        Amg::Vector3D dir = perigee->momentum().unit();

        const Amg::Transform3D* surfaceTransform = nullptr;
        const Amg::Transform3D* backupTransform = nullptr;
        std::unique_ptr<Amg::Transform3D> surfaceTransformToBeDeleted;
        double weightedDistanceSquared{0}, weightSquared{0};
        for (const Trk::TrackStateOnSurface* tsos : *states) {
            if (!tsos) continue;  // sanity check

            // require TrackParameters
            const Trk::TrackParameters* pars = tsos->trackParameters();
            if (!pars) continue;

            // check whether state is a measurement
            const Trk::MeasurementBase* meas = tsos->measurementOnTrack();
            if (!meas || tsos->type(Trk::TrackStateOnSurface::Outlier)) continue;
            rots->push_back(meas->clone());

            // only consider eta hits
            Identifier id = m_edmHelperSvc->getIdentifier(*meas);
            if (!id.is_valid() || m_idHelperSvc->measuresPhi(id)) continue;

            double distance = (pars->position() - perigee->position()).dot(dir);
            double weight = 1. / meas->localCovariance()(Trk::locX, Trk::locX);
            ATH_MSG_VERBOSE(" distance " << distance << " error " << Amg::error(meas->localCovariance(), Trk::locX) << " weight " << weight
                                         << " " << m_idHelperSvc->toString(id));
            weightedDistanceSquared += distance * weight;
            weightSquared += weight;
            if (m_idHelperSvc->isMdt(id)) {
                chIds.insert(m_idHelperSvc->chamberId(id));
                if (!surfaceTransform) {
                    const MdtDriftCircleOnTrack* mdt = dynamic_cast<const MdtDriftCircleOnTrack*>(meas);
                    if (mdt) {
                        // create new surface using AMDB reference frame
                        surfaceTransformToBeDeleted = std::make_unique<Amg::Transform3D>(mdt->detectorElement()->AmdbLRSToGlobalTransform().rotation());
                        surfaceTransformToBeDeleted->pretranslate(mdt->detectorElement()->center());
                        surfaceTransform = surfaceTransformToBeDeleted.get();                       
                    }
                }
            } else if ((m_idHelperSvc->isMM(id) || m_idHelperSvc->isCsc(id)) && !surfaceTransform) {
                surfaceTransform = &(meas)->associatedSurface().transform();
            } else if (!surfaceTransform && !backupTransform) {
                backupTransform = &(meas)->associatedSurface().transform();
            }
        }
        if (!surfaceTransform) surfaceTransform = backupTransform;
        // calculate distance new reference point, shift it 100 mm towards the start of the segment
        double refDistance = (weightSquared > 0 ? weightedDistanceSquared / weightSquared : 1) - 100;
        ATH_MSG_DEBUG(" weighted distance " << refDistance);

        const Amg::Vector3D refPos = perigee->position() + refDistance * dir;

        // find closest measured parameters
        double minDist = -1e6;
        const Trk::TrackParameters* closestPars = nullptr;
        for (const Trk::TrackStateOnSurface* tsos : *states) {
            if (!tsos) continue;  // sanity check

            // require TrackParameters
            const Trk::TrackParameters* pars = tsos->trackParameters();
            if (!pars || !pars->covariance()) continue;

            // look for the closest measured parameters to the reference point
            double distance = (pars->position() - refPos).dot(dir);
            if (distance < 0 && std::abs(distance) < std::abs(minDist)) {
                minDist = distance;
                closestPars = pars;
            }
        }

        if (!surfaceTransform) {
            ATH_MSG_WARNING(" failed to create a PlaneSurface for the track, cannot make segment!!! " << std::endl
                                                                                                      << m_printer->print(track)
                                                                                                      << std::endl
                                                                                                      << m_printer->printStations(track));
            return nullptr;
        }

        if (!closestPars) {
            closestPars = perigee;
            minDist = (perigee->position() - refPos).dot(dir);
        }

        Amg::Transform3D transform(surfaceTransform->rotation());
        transform.pretranslate(refPos);
        constexpr double surfDim = 500.;
        std::unique_ptr<Trk::PlaneSurface> surf = std::make_unique<Trk::PlaneSurface>(transform, surfDim, surfDim);
        std::unique_ptr<Trk::TrackParameters> exPars = m_propagator->propagate(ctx, *closestPars, *surf, minDist > 0 ? Trk::oppositeMomentum : Trk::alongMomentum, false,
                                              Trk::MagneticFieldProperties(Trk::NoField));
        if (!exPars || !exPars->covariance()) {
            ATH_MSG_VERBOSE("First trial reaching the surface failed. This is presumably due to a too large momentum. Let's try with a dummy 1 GeV momentum");
            std::unique_ptr<Trk::TrackParameters> cloned_pars {closestPars->clone()};
            constexpr double OneOverGeV = 1./ Gaudi::Units::GeV;
            cloned_pars->parameters()[Trk::qOverP] = cloned_pars->charge() *OneOverGeV;
            exPars = m_propagator->propagate(ctx, *cloned_pars, *surf,  Trk::anyDirection, false,
                                              Trk::MagneticFieldProperties(Trk::NoField));

            if (!exPars){
                ATH_MSG_DEBUG(" propagation failed!!! "<<*cloned_pars<<std::endl<<std::endl<<*surf);
                return nullptr;
            }
            exPars->parameters()[Trk::qOverP] = closestPars->parameters()[Trk::qOverP];
        }
        Amg::Vector2D locPos;
        if (!surf->globalToLocal(exPars->position(), exPars->momentum(), locPos)) {
            ATH_MSG_WARNING(" localToGlobal failed!!! ");            
            return nullptr;
        }
        Trk::LocalDirection locDir;
        surf->globalToLocalDirection(exPars->momentum(), locDir);

        // convert errors on global angles theta/phi to errors on local angles angleYZ/angleXZ
        Trk::JacobianPhiThetaLocalAngles globalToLocalMeasAnglesJacobian(exPars->parameters()[Trk::phi], exPars->parameters()[Trk::theta],
                                                                         exPars->associatedSurface().transform().rotation().inverse());

        // make the Jacobian to convert all in one go from global to local
        // so that the correlations are calculated correctly
        AmgSymMatrix(5) globalToLocalMeasJacobian;
        globalToLocalMeasJacobian.setZero();
        globalToLocalMeasJacobian(Trk::locX, Trk::locX) = 1.0;
        globalToLocalMeasJacobian(Trk::locY, Trk::locY) = 1.0;
        globalToLocalMeasJacobian(Trk::phi, Trk::phi) = globalToLocalMeasAnglesJacobian(0, 0);
        globalToLocalMeasJacobian(Trk::theta, Trk::theta) = globalToLocalMeasAnglesJacobian(1, 1);
        globalToLocalMeasJacobian(Trk::theta, Trk::phi) = globalToLocalMeasAnglesJacobian(0, 1);  // also fills (Trk::phi,Trk::theta)
        globalToLocalMeasJacobian(Trk::phi, Trk::theta) =  globalToLocalMeasJacobian(Trk::theta, Trk::phi);  // also fills (Trk::theta,Trk::phi)
        globalToLocalMeasJacobian(Trk::qOverP, Trk::qOverP) = 1.0;

        AmgSymMatrix(5) cov = exPars->covariance()->similarity(globalToLocalMeasJacobian);

        Trk::FitQuality* quality = nullptr;
        if (!chIds.empty()) {
            // calculate holes
            std::vector<Identifier> holes;
            for (const Identifier& chid : chIds) {
                std::vector<Identifier> holesChamber = calculateHoles(ctx, chid, *exPars, rots->stdcont());
                holes.insert(holes.end(), holesChamber.begin(), holesChamber.end());
            }
            quality = new MuonSegmentQuality(fq->chiSquared(), fq->numberDoF(), holes);

        } else {
            quality = new Trk::FitQuality(fq->chiSquared(), fq->numberDoF());
        }
        MuonSegment* seg = new MuonSegment(locPos, locDir, cov, surf.release(), rots.release(), quality);
        return seg;
    }

    std::vector<Identifier> MuonTrackToSegmentTool::calculateHoles(const EventContext& ctx, const Identifier& chid,
                                                                   const Trk::TrackParameters& pars,
                                                                   const MuonTrackToSegmentTool::MeasVec& measurements) const {
        SG::ReadCondHandle<Muon::MuonIntersectGeoData> InterSectSvc{m_chamberGeoKey,ctx};
        if (!InterSectSvc.isValid())   {
            ATH_MSG_ERROR("Failed to retrieve chamber intersection service");                            
        }

        const MuonStationIntersect intersect = InterSectSvc->tubesCrossedByTrack(chid, pars.position(), pars.momentum().unit());
        const MuonGM::MuonDetectorManager* MuonDetMgr = InterSectSvc->detMgr();

        // set to identify the hit on the segment
        std::set<Identifier> hitsOnSegment;
        for (const Trk::MeasurementBase* mit : measurements) {
            const MdtDriftCircleOnTrack* mdt = dynamic_cast<const MdtDriftCircleOnTrack*>(mit);
            if (mdt) hitsOnSegment.insert(mdt->identify());
        }

        // clear hole vector
        std::vector<Identifier> holes;
        for (unsigned int ii = 0; ii < intersect.tubeIntersects().size(); ++ii) {
            const MuonTubeIntersect& tint = intersect.tubeIntersects()[ii];

            // skip hole check if there is a hit in this tube
            if (hitsOnSegment.count(tint.tubeId)) continue;

            // if track goes through a tube which did not have a hit count as hole
            if (std::abs(tint.rIntersect) < MuonDetMgr->getMdtReadoutElement(tint.tubeId)->innerTubeRadius() && tint.xIntersect < -200.) {
                holes.push_back(tint.tubeId);
            }
        }
        return holes;
    }

}  // namespace Muon
