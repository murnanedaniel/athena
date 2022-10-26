/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonPatternSegmentMaker/MuonPatternCalibration.h"

#include <iostream>
#include <set>
#include <sstream>

#include "MuonPattern/MuonPatternChamberIntersect.h"
#include "MuonPrepRawData/MuonPrepDataContainer.h"
#include "MuonReadoutGeometry/CscReadoutElement.h"
#include "MuonReadoutGeometry/MdtReadoutElement.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/RpcReadoutElement.h"
#include "MuonReadoutGeometry/TgcReadoutElement.h"
#include "TrkToolInterfaces/IRIO_OnTrackCreator.h"

namespace Muon {

MuonPatternCalibration::MuonPatternCalibration(const std::string& t, const std::string& n, const IInterface* p)
    : AthAlgTool(t, n, p), m_keyRpc("RPC_Measurements"), m_keyTgc("TGC_Measurements")
{
    declareInterface<IMuonPatternCalibration>(this);

    declareProperty("DoMultiChamberAnalysis", m_doMultiAnalysis = true);
    declareProperty("DropDistance", m_dropDistance = 1500.);
    declareProperty("AngleCutPhi", m_phiAngleCut = 1e9);
    declareProperty("FullFinder", m_doFullFinder = true);
    declareProperty("RpcPrepDataContainer", m_keyRpc);
    declareProperty("TgcPrepDataContainer", m_keyTgc);
    declareProperty("DoSummary", m_doSummary = false);
    declareProperty("RecoverTriggerHits", m_recoverTriggerHits = true);
    declareProperty("RemoveDoubleMdtHits", m_removeDoubleMdtHits = true);
}

StatusCode
MuonPatternCalibration::initialize()
{
    ATH_MSG_VERBOSE("MuonPatternCalibration::Initializing");
    ATH_CHECK(m_mdtCreator.retrieve());
    ATH_CHECK(m_printer.retrieve());
    ATH_CHECK(m_idHelperSvc.retrieve());
    ATH_CHECK(m_clusterCreator.retrieve());
    ATH_CHECK(m_keyRpc.initialize());
    ATH_CHECK(m_keyTgc.initialize());
    return StatusCode::SUCCESS;
}


void
MuonPatternCalibration::calibrate(const MuonPatternCombination&           pattern,
                                  IMuonPatternCalibration::ROTsPerRegion& hitsPerRegion, const EventContext& ctx) const
{
    Muon::MuonPatternCalibration::RegionMap regionMap;
    bool                                    hasPhiMeasurements = checkForPhiMeasurements(pattern);
    createRegionMap(pattern, regionMap, hasPhiMeasurements, ctx);
    // calibrate hits
    calibrateRegionMap(regionMap, hitsPerRegion);
    return;
}


int
MuonPatternCalibration::getRegionId(const Identifier& id) const
{

    // simple division of MuonSpectrometer in regions using barrel/endcap seperation plus
    // inner/middle/outer seperation

    int stIndex = (int)m_idHelperSvc->stationIndex(id);
    int eta     = m_idHelperSvc->stationEta(id);

    int regionId = stIndex;
    if (eta < 0) regionId *= -1;

    return regionId;
}


bool
MuonPatternCalibration::checkForPhiMeasurements(const MuonPatternCombination& pat) const
{

    bool                                                     hasPhiMeasurements(false);
    std::vector<MuonPatternChamberIntersect>::const_iterator it = pat.chamberData().begin();
    for (; it != pat.chamberData().end(); ++it) {
        if (hasPhiMeasurements == true) {
            break;
        }
        std::vector<const Trk::PrepRawData*>::const_iterator pit = (*it).prepRawDataVec().begin();
        for (; pit != (*it).prepRawDataVec().end(); ++pit) {
            Identifier id = (*pit)->identify();
            if (m_idHelperSvc->isRpc(id) && m_idHelperSvc->rpcIdHelper().measuresPhi(id)) {
                hasPhiMeasurements = true;
                break;
            } else if (m_idHelperSvc->isTgc(id) && m_idHelperSvc->tgcIdHelper().isStrip(id)) {
                hasPhiMeasurements = true;
                break;
            } else if (m_idHelperSvc->isCsc(id) && m_idHelperSvc->cscIdHelper().measuresPhi(id)) {
                hasPhiMeasurements = true;
                break;
            }
        }
    }
    return hasPhiMeasurements;
}

void
MuonPatternCalibration::createRegionMap(const MuonPatternCombination& pat, RegionMap& regionMap,
                                        bool hasPhiMeasurements, const EventContext& ctx) const
{
    if (hasPhiMeasurements)
        ATH_MSG_DEBUG("pattern has phi measurements using extrapolation to determine second coordinate");
    else
        ATH_MSG_DEBUG("No phi measurements using center tubes");

    MuonPatternCalibration::Containers containers = retrieveTriggerHitContainers(ctx);

    std::vector<MuonPatternChamberIntersect>::const_iterator it = pat.chamberData().begin();
    for (; it != pat.chamberData().end(); ++it) {

        if (it->prepRawDataVec().empty()) continue;

        const Amg::Vector3D& patpose = (*it).intersectPosition();
        const Amg::Vector3D  patdire = (*it).intersectDirection().unit();

        /** Try to recover missing phi clusters:
            - loop over the clusters in the region and sort them by collection
            - count the number of eta and phi clusters per collection

        */
        std::map<int, EtaPhiHits> etaPhiHitsPerChamber;
        std::set<Identifier>      clusterIds;


        const Trk::PrepRawData* prd = (*it).prepRawDataVec().front();
        if (!prd) continue;
        Identifier id = prd->identify();

        // apply cut on the opening angle between pattern and chamber phi
        // do some magic to avoid problems at phi = 0 and 2*pi
        double phiStart  = patdire.phi();
        double chPhi     = prd->detectorElement()->center().phi();
        double pi        = M_PI;
        double phiRange  = 0.75 * pi;
        double phiRange2 = 0.25 * pi;
        double phiOffset = 0.;
        if (phiStart > phiRange || phiStart < -phiRange)
            phiOffset = 2 * pi;
        else if (phiStart > -phiRange2 && phiStart < phiRange2)
            phiOffset = pi;

        if (phiOffset > 1.5 * pi) {
            if (phiStart < 0) phiStart += phiOffset;
            if (chPhi < 0) chPhi += phiOffset;
        } else if (phiOffset > 0.) {
            phiStart += phiOffset;
            chPhi += phiOffset;
        }
        double dphi = fabs(phiStart - chPhi);

        if (dphi > m_phiAngleCut) {
            ATH_MSG_DEBUG("Large angular phi difference between pattern and chamber, phi pattern "
                          << patdire.phi() << "   phi chamber " << prd->detectorElement()->center().phi());
            continue;
        }

        // map to find duplicate hits in the chamber
        std::map<Identifier, const MdtPrepData*> idMdtMap;

        std::vector<const Trk::PrepRawData*>::const_iterator pit = (*it).prepRawDataVec().begin();
        for (; pit != (*it).prepRawDataVec().end(); ++pit) {

            if (!(*pit)) continue;

            const MdtPrepData* mdt = dynamic_cast<const MdtPrepData*>(*pit);
            if (mdt) {

                if (m_removeDoubleMdtHits) {
                    const MdtPrepData*& previousMdt = idMdtMap[mdt->identify()];
                    if (!previousMdt || previousMdt->tdc() > mdt->tdc())
                        previousMdt = mdt;
                    else
                        continue;
                }

                insertMdt(*mdt, regionMap, patpose, patdire, hasPhiMeasurements);
            } else {

                const MuonCluster* clus = dynamic_cast<const MuonCluster*>(*pit);
                if (!clus) continue;
                const Identifier& id = clus->identify();

                if (clusterIds.count(id)) continue;
                clusterIds.insert(id);

                if (m_recoverTriggerHits) {
                    bool measuresPhi = m_idHelperSvc->measuresPhi(id);
                    int  colHash     = clus->collectionHash();

                    EtaPhiHits& hitsPerChamber = etaPhiHitsPerChamber[colHash];
                    if (measuresPhi)
                        ++hitsPerChamber.nphi;
                    else
                        ++hitsPerChamber.neta;
                }
                insertCluster(*clus, regionMap, patpose, patdire, hasPhiMeasurements);
            }
        }
        if (!etaPhiHitsPerChamber.empty()) {
            std::map<int, EtaPhiHits>::iterator chit     = etaPhiHitsPerChamber.begin();
            std::map<int, EtaPhiHits>::iterator chit_end = etaPhiHitsPerChamber.end();
            for (; chit != chit_end; ++chit) {
                EtaPhiHits& hits = chit->second;
                if ((hits.neta > 0 && hits.nphi == 0) || (hits.nphi > 0 && hits.neta == 0)) {
                    if (m_idHelperSvc->isRpc(id) && containers.m_rpcPrdContainer) {

                        auto pos = containers.m_rpcPrdContainer->indexFindPtr(chit->first);
                        if (pos == nullptr)
                            ATH_MSG_DEBUG("RpcPrepDataCollection not found in container!!"<< m_keyRpc);
                        else {
                            RpcPrepDataCollection::const_iterator rpcit     = pos->begin();
                            RpcPrepDataCollection::const_iterator rpcit_end = pos->end();
                            for (; rpcit != rpcit_end; ++rpcit) {
                                if (clusterIds.count((*rpcit)->identify())) continue;
                                const MuonCluster* clus = dynamic_cast<const MuonCluster*>(*rpcit);
                                insertCluster(*clus, regionMap, patpose, patdire, hasPhiMeasurements);
                                clusterIds.insert(clus->identify());
                            }
                        }
                    } else if (m_idHelperSvc->isTgc(id) && containers.m_tgcPrdContainer) {
                        auto pos = containers.m_tgcPrdContainer->indexFindPtr(chit->first);
                        if (pos == nullptr)
                            ATH_MSG_DEBUG("TgcPrepDataCollection not found in container!! "<< m_keyTgc);
                        else {
                            TgcPrepDataCollection::const_iterator tgcit     = pos->begin();
                            TgcPrepDataCollection::const_iterator tgcit_end = pos->end();
                            for (; tgcit != tgcit_end; ++tgcit) {
                                if (clusterIds.count((*tgcit)->identify())) continue;
                                const MuonCluster* clus = dynamic_cast<const MuonCluster*>(*tgcit);
                                insertCluster(*clus, regionMap, patpose, patdire, hasPhiMeasurements);
                                clusterIds.insert(clus->identify());
                            }
                        }
                    }
                }
            }
        }
    }
}

void
MuonPatternCalibration::insertCluster(const MuonCluster& clus, RegionMap& regionMap, const Amg::Vector3D& patpose,
                                      const Amg::Vector3D& patdire, bool hasPhiMeasurements) const
{

    const Identifier& id = clus.identify();
    // check whether we are measuring phi or eta
    bool measuresPhi = m_idHelperSvc->measuresPhi(id);

    Amg::Vector3D globalpos = clus.globalPosition();
    Amg::Vector3D intersect;

    if (hasPhiMeasurements) {
        // if there is a phi measurement in the pattern use the global direction to calculate the intersect with
        // measurement plane and use the intersect to calculate the position along the strip

        // calculate intersect pattern measurement plane
        const Trk::Surface& surf         = clus.detectorElement()->surface(id);
        Amg::Vector3D       planepostion = surf.center();
        Amg::Vector3D       planenormal  = surf.normal();
        double              denom        = patdire.dot(planenormal);
        double              u            = (planenormal.dot(planepostion - patpose)) / (denom);
        Amg::Vector3D       piOnPlane    = (patpose + u * patdire);

        // transform to local plane coordiantes
        const Amg::Transform3D gToLocal = clus.detectorElement()->surface().transform().inverse();
        Amg::Vector3D          ilpos    = gToLocal * piOnPlane;
        Amg::Vector3D          glpos    = gToLocal * globalpos;

        // strip length
        double striplen(0.);

        // projective strips
        bool hasPointingPhiStrips = false;

        // detector specific stuff
        const RpcPrepData* rpc = dynamic_cast<const RpcPrepData*>(&clus);
        if (rpc) {
            striplen = rpc->detectorElement()->StripLength(measuresPhi);
        } else {
            const TgcPrepData* tgc = dynamic_cast<const TgcPrepData*>(&clus);
            if (!tgc) return;

            int gasGap = m_idHelperSvc->tgcIdHelper().gasGap(id);
            if (measuresPhi) {
                hasPointingPhiStrips = true;
                striplen             = tgc->detectorElement()->StripLength(gasGap);
            } else {
                int wire = m_idHelperSvc->tgcIdHelper().channel(id);
                striplen = tgc->detectorElement()->WireLength(gasGap, wire);
            }
        }

        // set the position along the strip
        if (!measuresPhi) {
            glpos[0] = ilpos.x();
        } else {
            if (hasPointingPhiStrips) {
                // do some special for tgcs
                glpos[1] = ilpos.y();
            } else {
                glpos[1] = ilpos.y();
            }
        }

        // transform back to global coordinates
        intersect         = gToLocal.inverse() * glpos;
        Amg::Vector3D dif = globalpos - intersect;
        if ((intersect - piOnPlane).mag() > m_dropDistance || dif.mag() > 0.5 * striplen + m_dropDistance) {

            ATH_MSG_VERBOSE(">>>> extrapolated position far outside volume, dropping hit "
                            << m_idHelperSvc->toString(id) << ". dist along strip " << dif.mag() << " 1/2 strip len "
                            << 0.5 * striplen << " dist measurement plane " << (intersect - piOnPlane).mag());
            return;
        }
        if (dif.mag() > 0.5 * striplen) {
            Amg::Vector3D newpos = globalpos - dif * (0.5 * striplen / dif.mag());

            ATH_MSG_VERBOSE(">>>> extrapolated position outside volume, shifting position "
                            << m_idHelperSvc->toString(id) << ". position along strip " << dif.mag() << " 1/2 tube len "
                            << 0.5 * striplen << " dist To strip " << (intersect - piOnPlane).mag()
                            << ". dist to newpos " << (newpos - globalpos).mag() << " pos " << newpos);

            intersect = newpos;
        }
    } else {
        // no phi measurements, use strip center
        intersect = globalpos;
    }

    // enter hit in map
    int regionId = getRegionId(id);

    RegionMapIt mip = regionMap.find(regionId);
    if (mip == regionMap.end()) {
        Region region;
        region.triggerPrds.push_back(std::make_pair(intersect, &clus));
        region.regionDir = patdire;
        region.regionPos = patpose;
        regionMap.insert(std::make_pair(regionId, region));
    } else {
        Region& region = mip->second;
        region.triggerPrds.push_back(std::make_pair(intersect, &clus));
    }
}


void
MuonPatternCalibration::insertMdt(const MdtPrepData& mdt, RegionMap& regionMap, const Amg::Vector3D& patpose,
                                  const Amg::Vector3D& patdire, bool hasPhiMeasurements) const
{

    Amg::Vector3D     intersect;
    const Identifier& id = mdt.identify();

    const MuonGM::MdtReadoutElement* detEl   = mdt.detectorElement();
    const Amg::Vector3D&             tubePos = mdt.globalPosition();

    if (hasPhiMeasurements) {
        // if there is a phi measurement in the pattern use the global direction to calculate the intersect with the
        // tube use the intersect to calculate the second coordinate

        const Amg::Transform3D amdbToGlobal = detEl->AmdbLRSToGlobalTransform();


        // calculate intersect pattern measurement plane
        const Amg::Vector3D& planepostion = tubePos;

        // always project on plane with normal in radial direction
        Amg::Vector3D planenormal = m_idHelperSvc->mdtIdHelper().isBarrel(id)
                                        ? amdbToGlobal.linear() * Amg::Vector3D(0., 0., 1.)
                                        : amdbToGlobal.linear() * Amg::Vector3D(0., 1., 0.);

        double        denom     = patdire.dot(planenormal);
        double        u         = (planenormal.dot(planepostion - patpose)) / (denom);
        Amg::Vector3D piOnPlane = (patpose + u * patdire);

        Amg::Vector3D lpiOnPlane = amdbToGlobal.inverse() * piOnPlane;
        Amg::Vector3D ltubePos   = amdbToGlobal.inverse() * tubePos;

        intersect = amdbToGlobal * Amg::Vector3D(lpiOnPlane.x(), ltubePos.y(), ltubePos.z());

        Amg::Vector3D dif     = tubePos - intersect;
        double        tubelen = detEl->getActiveTubeLength(m_idHelperSvc->mdtIdHelper().tubeLayer(id),
                                                    m_idHelperSvc->mdtIdHelper().tube(id));

        if (dif.mag() > 0.5 * tubelen) {
            Amg::Vector3D newpos = tubePos - dif * (0.5 * tubelen / dif.mag());

            ATH_MSG_VERBOSE(">>>> extrapolated position outside volume, shifting position "
                            << m_idHelperSvc->toString(id) << ". position along strip " << dif.mag() << " 1/2 tube len "
                            << 0.5 * tubelen << " dist To Wire " << (piOnPlane - intersect).mag() << ". dist to newpos "
                            << (newpos - tubePos).mag() << " pos " << newpos);

            intersect = newpos;
        }
    } else {
        // not phi measurement, use tube center
        intersect = tubePos;
    }

    // enter hit in map
    Identifier elId = m_idHelperSvc->mdtIdHelper().elementID(id);

    MuonStationIndex::ChIndex chIndex = m_idHelperSvc->chamberIndex(elId);
    int                       chFlag  = elId.get_identifier32().get_compact();
    if (m_doMultiAnalysis) {
        if (m_idHelperSvc->isSmallChamber(id)) {
            ATH_MSG_VERBOSE("Small chamber " << m_idHelperSvc->toString(elId));
            chFlag = 0;
            if (chIndex == MuonStationIndex::BIS) {
                int eta = m_idHelperSvc->stationEta(elId);
                if (abs(eta) == 8) {
                    ATH_MSG_VERBOSE("BIS8 chamber " << m_idHelperSvc->toString(elId));
                    chFlag = 3;
                }
            }
        } else {
            ATH_MSG_VERBOSE("Large chamber " << m_idHelperSvc->toString(elId));
            chFlag = 1;
            if (chIndex == MuonStationIndex::BIL) {
                std::string stName = m_idHelperSvc->chamberNameString(id);
                if (stName[2] == 'R') {
                    ATH_MSG_VERBOSE("BIR chamber " << m_idHelperSvc->toString(elId));
                    chFlag = 2;
                }
            }
        }
        int phi = m_idHelperSvc->mdtIdHelper().stationPhi(id);

        chFlag += 10 * phi;
    }
    // use center tube for region assignment
    int regionId = getRegionId(id);

    RegionMapIt mip = regionMap.find(regionId);
    if (mip == regionMap.end()) {
        Region region;
        region.mdtPrdsPerChamber[chFlag].push_back(std::make_pair(intersect, &mdt));
        region.regionPos = patpose;
        region.regionDir = patdire;
        regionMap.insert(std::make_pair(regionId, region));
    } else {
        Region& region = mip->second;
        region.mdtPrdsPerChamber[chFlag].push_back(std::make_pair(intersect, &mdt));
    }
}

void
MuonPatternCalibration::clearRotsPerRegion(IMuonPatternCalibration::ROTsPerRegion& hitsPerRegion) const
{
    // loop over regions, delete all ROTS in map
    ROTsPerRegionIt rit     = hitsPerRegion.begin();
    ROTsPerRegionIt rit_end = hitsPerRegion.end();
    for (; rit != rit_end; ++rit) {

        // loop over mdt hits vectors in region
        MdtVecIt mdtvit     = rit->mdts.begin();
        MdtVecIt mdtvit_end = rit->mdts.end();
        for (; mdtvit != mdtvit_end; ++mdtvit) {

            // loop over mdt hits
            MdtIt mdtit     = mdtvit->begin();
            MdtIt mdtit_end = mdtvit->end();
            for (; mdtit != mdtit_end; ++mdtit) {
                delete *mdtit;
            }
        }

        // loop over clusters in region
        ClusterIt clit     = rit->clusters.begin();
        ClusterIt clit_end = rit->clusters.end();
        for (; clit != clit_end; ++clit) {
            delete *clit;
        }
    }
}


void
MuonPatternCalibration::printRegionMap(const RegionMap& regionMap) const
{

    RegionMapCit mit     = regionMap.begin();
    RegionMapCit mit_end = regionMap.end();
    ATH_MSG_INFO("Summarizing input");

    for (; mit != mit_end; ++mit) {
        ATH_MSG_INFO("new region " << mit->first << " trigger " << mit->second.triggerPrds.size() << " mdt ch "
                                   << mit->second.mdtPrdsPerChamber.size());

        if (!mit->second.triggerPrds.empty()) ATH_MSG_INFO("trigger hits " << mit->second.triggerPrds.size());


        ISPrdIt pit     = mit->second.triggerPrds.begin();
        ISPrdIt pit_end = mit->second.triggerPrds.end();
        for (; pit != pit_end; ++pit) {
            ATH_MSG_INFO("  " << m_printer->print(*(pit->second)));
        }

        RegionIdMapIt idit     = mit->second.mdtPrdsPerChamber.begin();
        RegionIdMapIt idit_end = mit->second.mdtPrdsPerChamber.end();
        for (; idit != idit_end; ++idit) {
            ATH_MSG_INFO("new MDT chamber with " << idit->second.size() << " hits");

            ISPrdMdtIt mdtit     = idit->second.begin();
            ISPrdMdtIt mdtit_end = idit->second.end();
            for (; mdtit != mdtit_end; ++mdtit) {
                const MdtPrepData* prd = mdtit->second;
                ATH_MSG_INFO("  " << m_printer->print(*prd));
            }
        }
    }
}

void
MuonPatternCalibration::calibrateRegionMap(const RegionMap&                        regionMap,
                                           IMuonPatternCalibration::ROTsPerRegion& hitsPerRegion) const
{

    RegionMapCit mit     = regionMap.begin();
    RegionMapCit mit_end = regionMap.end();

    for (; mit != mit_end; ++mit) {


        ROTRegion rotRegion;
        rotRegion.regionId  = mit->first;
        rotRegion.regionPos = mit->second.regionPos;
        rotRegion.regionDir = mit->second.regionDir;


        ISPrdIt pit     = mit->second.triggerPrds.begin();
        ISPrdIt pit_end = mit->second.triggerPrds.end();
        for (; pit != pit_end; ++pit) {

            if (!pit->second) continue;

            const MuonClusterOnTrack* cluster = m_clusterCreator->createRIO_OnTrack(*(pit->second), pit->first);
            if (!cluster) continue;

            rotRegion.clusters.push_back(cluster);
        }

        RegionIdMapIt idit     = mit->second.mdtPrdsPerChamber.begin();
        RegionIdMapIt idit_end = mit->second.mdtPrdsPerChamber.end();
        for (; idit != idit_end; ++idit) {

            ISPrdMdtIt mdtit     = idit->second.begin();
            ISPrdMdtIt mdtit_end = idit->second.end();

            MdtVec mdtROTs;
            for (; mdtit != mdtit_end; ++mdtit) {
                const MdtPrepData* prd = mdtit->second;
                if (!prd) continue;

                const MdtDriftCircleOnTrack* mdt = m_mdtCreator->createRIO_OnTrack(*prd, mdtit->first);

                if (!mdt) {
                    ATH_MSG_VERBOSE("Failed to calibrate " << m_idHelperSvc->toString(prd->identify()));
                    continue;
                }
                mdtROTs.push_back(mdt);
            }
            if (!mdtROTs.empty()) rotRegion.mdts.push_back(mdtROTs);
        }
        hitsPerRegion.push_back(rotRegion);
    }
}

MuonPatternCalibration::Containers
MuonPatternCalibration::retrieveTriggerHitContainers(const EventContext& ctx) const
{
    MuonPatternCalibration::Containers containers;
    SG::ReadHandle<Muon::RpcPrepDataContainer> RpcCont(m_keyRpc, ctx);
    if (!RpcCont.isValid()) {
        ATH_MSG_DEBUG(" Failed to retrieve RpcPrepDataContainer, will not recover rpc trigger hits ");
    } else{
        containers.m_rpcPrdContainer = RpcCont.cptr();
    }
    SG::ReadHandle<Muon::TgcPrepDataContainer> TgcCont(m_keyTgc, ctx);
    if (!TgcCont.isValid()) {
        ATH_MSG_DEBUG(" Failed to retrieve TgcPrepDataContainer, will not recover tgc trigger hits ");
    } else{
        containers.m_tgcPrdContainer = TgcCont.cptr();
    }

    return containers;
}

}  // namespace Muon
