/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// Muon
#include "MuonTrackingGeometry/MuonInertMaterialBuilder.h"

#include "MuonTrackingGeometry/MuonStationTypeBuilder.h"
// MuonSpectrometer include
#include "MuonReadoutGeometry/MuonDetectorManager.h"
// Amg
#include "GeoPrimitives/CLHEPtoEigenConverter.h"
#include "GeoPrimitives/GeoPrimitives.h"
// Trk
#include "TrkDetDescrInterfaces/IDetachedTrackingVolumeBuilder.h"
#include "TrkDetDescrInterfaces/ILayerArrayCreator.h"
#include "TrkDetDescrInterfaces/ILayerBuilder.h"
#include "TrkDetDescrInterfaces/ITrackingVolumeArrayCreator.h"
#include "TrkDetDescrUtils/GeometryStatics.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkGeometry/CylinderLayer.h"
#include "TrkGeometry/DiscLayer.h"
#include "TrkGeometry/PlaneLayer.h"
#include "TrkGeometry/SubtractedCylinderLayer.h"
#include "TrkGeometry/SubtractedPlaneLayer.h"
#include "TrkGeometry/TrackingGeometry.h"
#include "TrkGeometry/TrackingVolume.h"
#include "TrkSurfaces/CylinderSurface.h"
#include "TrkSurfaces/DiscBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
#include "TrkSurfaces/TrapezoidBounds.h"
#include "TrkVolumes/BoundarySurface.h"
#include "TrkVolumes/CombinedVolumeBounds.h"
#include "TrkVolumes/CuboidVolumeBounds.h"
#include "TrkVolumes/PrismVolumeBounds.h"
#include "TrkVolumes/SimplePolygonBrepVolumeBounds.h"
#include "TrkVolumes/SubtractedVolumeBounds.h"
#include "TrkVolumes/TrapezoidVolumeBounds.h"
#include "TrkVolumes/VolumeExcluder.h"

// STD
#include <fstream>
#include <iostream>
#include <map>

#include "GeoModelKernel/GeoCons.h"
#include "GeoModelKernel/GeoShape.h"
#include "GeoModelKernel/GeoShapeShift.h"
#include "GeoModelKernel/GeoTube.h"
#include "GeoModelKernel/GeoTubs.h"
// mw
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoPara.h"
#include "GeoModelKernel/GeoPgon.h"
#include "GeoModelKernel/GeoShapeIntersection.h"
#include "GeoModelKernel/GeoShapeSubtraction.h"
#include "GeoModelKernel/GeoShapeUnion.h"
#include "GeoModelKernel/GeoTrap.h"
#include "GeoModelKernel/GeoTrd.h"
#include "GeoModelKernel/GeoVolumeCursor.h"
#include "TrkSurfaces/EllipseBounds.h"

// constructor
Muon::MuonInertMaterialBuilder::MuonInertMaterialBuilder(const std::string& t, const std::string& n, const IInterface* p) :
    AthAlgTool(t, n, p), Trk::TrackingVolumeManipulator() {
    declareInterface<Trk::IDetachedTrackingVolumeBuilder>(this);
}

// Athena standard methods
// initialize
StatusCode Muon::MuonInertMaterialBuilder::initialize() {
    StatusCode ds = detStore()->retrieve(m_muonMgr);

    if (ds.isFailure()) { ATH_MSG_ERROR("Could not get MuonDetectorManager, no layers for muons will be built. "); }

    ATH_MSG_INFO(m_muonMgr->geometryVersion());

    // if no muon materials are declared, take default ones

    m_muonMaterial = Trk::Material(10e10, 10e10, 0., 0., 0.);  // default material properties

    // mw
    m_materialConverter = std::make_unique<Trk::GeoMaterialConverter>();
    if (!m_materialConverter) {
        ATH_MSG_FATAL("Could not create material converter in " << name() << " initialize()");
        return StatusCode::FAILURE;
    }
    m_geoShapeConverter = std::make_unique<Trk::GeoShapeConverter>();
    if (!m_geoShapeConverter) {
        ATH_MSG_FATAL("Could not create shape converter in " << name() << " initialize()");
        return StatusCode::FAILURE;
    }

    // random number generator
    if (m_rndmGenSvc.retrieve().isFailure()) {
        ATH_MSG_FATAL("Could not retrieve " << m_rndmGenSvc);
        return StatusCode::FAILURE;
    }
    m_flatDist = std::make_unique<Rndm::Numbers>(&(*m_rndmGenSvc), Rndm::Flat(0., 1.));
    if (!m_flatDist) {
        ATH_MSG_FATAL("Could not create flat distribution random generator in " << name() << " initialize()");
        return StatusCode::FAILURE;
    }

    if (m_simplifyToLayers) { ATH_MSG_INFO(name() << " option Simplify(Muon)GeometryToLayers no longer maintained "); }

    m_constituents.clear();

    ATH_MSG_INFO(name() << " initialize() successful");

    return StatusCode::SUCCESS;
}

const std::vector<const Trk::DetachedTrackingVolume*>* Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumes(bool blend) {
    // split output into objects to be kept and objects which may be released from memory (blended)
    std::pair<std::vector<const Trk::DetachedTrackingVolume*>, std::vector<const Trk::DetachedTrackingVolume*> > mInert;

    // retrieve muon station prototypes from GeoModel
    const std::vector<std::pair<const Trk::DetachedTrackingVolume*, std::vector<Amg::Transform3D> > >* msTypes =
        buildDetachedTrackingVolumeTypes(blend);
    ATH_MSG_INFO(name() << " obtained " << msTypes->size() << " prototypes");

    std::vector<std::pair<const Trk::DetachedTrackingVolume*, std::vector<Amg::Transform3D> > >::const_iterator msTypeIter =
        msTypes->begin();

    for (; msTypeIter != msTypes->end(); ++msTypeIter) {
        std::string msTypeName = (*msTypeIter).first->name();
        bool perm = true;
        if (blend) {
            // decide if object suitable for blending; does not concern shields
            double protMass = 0.;
            for (unsigned int ic = 0; ic < (*msTypeIter).first->constituents()->size(); ic++) {
                protMass += calculateVolume((*((*msTypeIter).first->constituents()))[ic].first.get()) *
                            (*((*msTypeIter).first->constituents()))[ic].second;
            }
            perm = msTypeName.compare(0, 1, "J") != 0 && m_blendLimit > 0 && protMass > m_blendLimit;
        }
        if (perm) msTypeName +=  "PERM";
        //
        const Trk::DetachedTrackingVolume* msTV = (*msTypeIter).first;
        for (unsigned int it = 0; it < (*msTypeIter).second.size(); it++) {
            Amg::Transform3D combTr((*msTypeIter).second[it]);
            const Trk::DetachedTrackingVolume* newStat = msTV->clone(msTypeName, combTr);
            if (perm)
                mInert.first.push_back(newStat);
            else
                mInert.second.push_back(newStat);
        }
    }

    // clean up prototypes
    for (unsigned int it = 0; it < msTypes->size(); it++) delete (*msTypes)[it].first;
    delete msTypes;

    // merge
    const std::vector<const Trk::DetachedTrackingVolume*>* muonObjects = nullptr;
    if (mInert.first.empty())
        muonObjects = new std::vector<const Trk::DetachedTrackingVolume*>(mInert.second);
    else {
        for (unsigned int i = 0; i < mInert.second.size(); i++) mInert.first.push_back(mInert.second[i]);
        muonObjects = new std::vector<const Trk::DetachedTrackingVolume*>(mInert.first);
    }

    ATH_MSG_INFO(name() << " returns  " << (*muonObjects).size() << " objects (detached volumes)");

    return muonObjects;
}

const std::vector<std::pair<const Trk::DetachedTrackingVolume*, std::vector<Amg::Transform3D> > >*
Muon::MuonInertMaterialBuilder::buildDetachedTrackingVolumeTypes(bool blend) {
    ATH_MSG_INFO(name() << " building muon object types");
    std::vector<std::pair<const Trk::DetachedTrackingVolume*, std::vector<Amg::Transform3D> > > objs;

    std::vector<std::string> objName;
    if (!m_muonMgr) {
        ATH_MSG_FATAL("Passive material cannot be built with a detector manager");
        return nullptr;
    }
    // link to top tree
    const GeoVPhysVol* top = &(*(m_muonMgr->getTreeTop(0)));
    GeoVolumeCursor vol(top);

    while (!vol.atEnd()) {
        const GeoVPhysVol* cv = &(*(vol.getVolume()));
        const GeoLogVol* clv = cv->getLogVol();
        const std::string vnamestr = clv->getName();
        std::string_view vname(vnamestr);
        // NSW station
        if (vname.size() > 7 && vname.substr(vname.size() - 7, 7) == "Station") {  // do nothing, active station
        } else {
            // std::cout << " INERT muon object found:" << vname << std::endl;

            bool accepted = true;
            if (vname.substr(0, 3) == "BAR" || vname.substr(0, 2) == "BT" || vname.substr(0, 6) == "EdgeBT" ||
                vname.substr(0, 6) == "HeadBT")
                accepted = m_buildBT;
            else if (vname.substr(0, 3) == "ECT")
                accepted = m_buildECT;
            else if (vname.substr(0, 4) == "Feet" || (vname.size() > 7 && (vname.substr(3, 4) == "Feet" || vname.substr(4, 4) == "Feet")))
                accepted = m_buildFeets;
            else if (vname.substr(0, 4) == "Rail")
                accepted = m_buildRails > 0;
            else if (vname.substr(0, 1) == "J")
                accepted = m_buildShields > 0;
            // NSW build inertmaterial for spacer frame, aluminium HUB, NJD disk and A plate
            else if (vname.substr(0, 3) == "NSW" && vname.substr(1, 6) == "Spacer")
                accepted = m_buildNSWInert;
            else if (vname.substr(0, 3) == "NSW" && vname.substr(1, 2) == "Al")
                accepted = m_buildNSWInert;
            else if (vname.substr(0, 3) == "NJD")
                accepted = m_buildNSWInert;
            else if (vname.substr(0, 1) == "A" && vname.substr(1, 5) == "Plate")
                accepted = m_buildNSWInert;
            // strange NSW will be anyway build
            else if (vname.substr(0, 1) != "J")
                accepted = m_buildSupports > 0;

            // if ( vname=="EdgeBTVoussoir" && accepted && m_simplify ) accepted = false;

            if (accepted) ATH_MSG_VERBOSE(name() << " INERT muon object found:" << vname);
            if (accepted) ATH_MSG_VERBOSE(" INERT muon object found and accepted :" << vname);
            if (!accepted) ATH_MSG_VERBOSE(" INERT muon object found and rejected :" << vname);

            if (!accepted) {
                vol.next();
                continue;
            }

            // update to accomodate AGDD structures

            if (msg().level() == MSG::VERBOSE) printInfo(cv);

            std::vector<const GeoShape*> input_shapes;
            std::vector<std::pair<const GeoLogVol*, std::vector<Amg::Transform3D> > > vols;

            bool simpleTree = false;
            if (!cv->getNChildVols()) {
                std::vector<Amg::Transform3D> volTr;
                volTr.push_back(vol.getTransform());
                vols.emplace_back(clv, std::move(volTr));
                simpleTree = true;
            } else {
                getObjsForTranslation(cv, Trk::s_idTransform, vols);
            }
            input_shapes.resize(vols.size());
            for (unsigned int i = 0; i < vols.size(); i++) input_shapes[i] = vols[i].first->getShape();

            // reserve m_constituents to avoid having DetachedTrackingVolume.constituents() point to a Vector that was reallocated
            // meanwhile
            m_constituents.reserve(vols.size() + 3);

            for (unsigned int ish = 0; ish < vols.size(); ish++) {
                std::string protoName = std::string(vname);
                if (!simpleTree) protoName += (vols[ish].first->getName());
                std::string pName = vols[ish].first->getName();
                ATH_MSG_VERBOSE(" check in pName " << pName << ", made of " << vols[ish].first->getMaterial()->getName() << " x0 "
                                                   << vols[ish].first->getMaterial()->getRadLength() << ","
                                                   << vols[ish].first->getShape()->type());
                // if(pName.substr(0,6)=="sTGC_1") continue;
                // if(pName.substr(0,4)=="MM_1") continue;

                bool found = false;
                for (unsigned int ip = 0; ip < objs.size(); ip++) {
                    if (protoName == objs[ip].first->name()) {
                        found = true;
                        if (simpleTree) objs[ip].second.push_back(vol.getTransform());
                        // else objs[ip].second.insert(objs[ip].second.end(),vols[ish].second.begin(),vols[ish].second.end());
                    }
                }
                if (found) continue;
                // m_geoShapeConverter->decodeShape(input_shapes[ish]);
                Amg::Transform3D ident;
                ident.setIdentity();
                if (msg().level() == MSG::VERBOSE) {
                    const Trk::Volume* trTest = m_geoShapeConverter->translateGeoShape(input_shapes[ish], &(vols[ish].second[0]));
                    delete trTest;
                }
                const Trk::Volume* trObject = m_geoShapeConverter->translateGeoShape(input_shapes[ish], &ident);

                if (trObject) {
                    Trk::Material mat = m_materialConverter->convert(vols[ish].first->getMaterial());
                    Trk::TrackingVolume* newType = new Trk::TrackingVolume(*trObject, mat, nullptr, nullptr, protoName);
                    Trk::TrackingVolume* simType = simplifyShape(newType, blend);
                    Trk::DetachedTrackingVolume* typeStat = new Trk::DetachedTrackingVolume(protoName, simType);
                    if (blend) typeStat->saveConstituents(&(m_constituents.back()));
                    objs.emplace_back(typeStat, vols[ish].second);
                    delete trObject;

                } else {
                    ATH_MSG_WARNING(name() << " volume not translated: " << vname);
                }
            }  // end new object
        }
        vol.next();
    }

    // add some extra material to the endcap     // ST this should be phased out eventually ?
    if (m_extraMaterial) {
        // scaling
        float scmat1 = m_extraX0 * m_extraFraction / 10. * 88.93;
        float scmat2 = m_extraX0 * (1. - m_extraFraction) / 10. * 88.93;
        // used : z = 14; A=28 ; rho = 2.33 g/cm^3, X0 = 93.7 mmm, l0 = 465.2 mm (Silicium)
        Trk::Material mat1(93.7 / scmat1, 465.2 / scmat1, scmat1 * 14, scmat1 * 28, 0.0023, 0.);
        Trk::Material mat2(93.7 / scmat2, 465.2 / scmat2, scmat2 * 14, scmat2 * 28, 0.0023, 0.);
        Trk::LayerArray* dummyLayers = nullptr;
        Trk::TrackingVolumeArray* dummyVolumes = nullptr;
        Trk::VolumeBounds* extraBounds1 = new Trk::CylinderVolumeBounds(850., 13000., 5.);
        Trk::TrackingVolume* mextra1 = new Trk::TrackingVolume(nullptr, extraBounds1, mat1, dummyLayers, dummyVolumes, "extraMat1");
        Trk::TrackingVolume* simType1 = simplifyShape(mextra1, blend);
        Trk::DetachedTrackingVolume* eVol1 = new Trk::DetachedTrackingVolume("extraTGCmat1", simType1);
        if (blend) eVol1->saveConstituents(&(m_constituents.back()));
        Trk::VolumeBounds* extraBounds2 = new Trk::CylinderVolumeBounds(850., 13000., 5.);
        Trk::TrackingVolume* mextra2 = new Trk::TrackingVolume(nullptr, extraBounds2, mat2, dummyLayers, dummyVolumes, "extraMat2");
        Trk::TrackingVolume* simType2 = simplifyShape(mextra2, blend);
        Trk::DetachedTrackingVolume* eVol2 = new Trk::DetachedTrackingVolume("extraTGCmat2", simType2);
        if (blend) eVol2->saveConstituents(&(m_constituents.back()));
        std::vector<Amg::Transform3D> pos1;
        pos1.emplace_back(Amg::Translation3D(0., 0., m_extraPos1));
        pos1.emplace_back(Amg::Translation3D(0., 0., -m_extraPos1));
        std::vector<Amg::Transform3D> pos2;
        pos2.emplace_back(Amg::Translation3D(0., 0., m_extraPos2));
        pos2.emplace_back(Amg::Translation3D(0., 0., -m_extraPos2));
        objs.emplace_back(eVol1, pos1);
        objs.emplace_back(eVol2, pos2);
    }
    //

    const std::vector<std::pair<const Trk::DetachedTrackingVolume*, std::vector<Amg::Transform3D> > >* mObjects =
        new std::vector<std::pair<const Trk::DetachedTrackingVolume*, std::vector<Amg::Transform3D> > >(objs);

    int count = 0;
    for (unsigned int i = 0; i < mObjects->size(); i++) count += (*mObjects)[i].second.size();

    ATH_MSG_INFO(name() << " returns " << mObjects->size() << " prototypes, to be cloned into " << count << " objects");

    return mObjects;
}

// finalize
StatusCode Muon::MuonInertMaterialBuilder::finalize() {
    ATH_MSG_INFO(name() << " finalize() successful");
    return StatusCode::SUCCESS;
}
//

void Muon::MuonInertMaterialBuilder::printInfo(const GeoVPhysVol* pv) const {
    const GeoLogVol* lv = pv->getLogVol();
    ATH_MSG_VERBOSE("New Muon Inert Object:" << lv->getName() << ", made of " << lv->getMaterial()->getName() << " x0 "
                                             << lv->getMaterial()->getRadLength() << "," << lv->getShape()->type());
    m_geoShapeConverter->decodeShape(lv->getShape());
    printChildren(pv);
}

void Muon::MuonInertMaterialBuilder::printChildren(const GeoVPhysVol* pv) const {
    // subcomponents
    unsigned int nc = pv->getNChildVols();
    for (unsigned int ic = 0; ic < nc; ic++) {
        Amg::Transform3D transf = pv->getXToChildVol(ic);

        //
        // std::cout << " dumping transform to subcomponent" << std::endl;
        // std::cout << transf[0][0]<<"," <<transf[0][1]<<"," <<transf[0][2]<<","<<transf[0][3] << std::endl;
        // std::cout << transf[1][0]<<"," <<transf[1][1]<<"," <<transf[1][2]<<","<<transf[1][3] << std::endl;
        // std::cout << transf[2][0]<<"," <<transf[2][1]<<"," <<transf[2][2]<<","<<transf[2][3] << std::endl;
        //
        const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
        const GeoLogVol* clv = cv->getLogVol();
        ATH_MSG_VERBOSE("  ");
        ATH_MSG_VERBOSE("subcomponent:" << ic << ":" << clv->getName() << ", made of " << clv->getMaterial()->getName() << " x0 "
                                        << clv->getMaterial()->getRadLength() << " , " << clv->getShape()->type() << ","
                                        << transf.translation().x() << " " << transf.translation().y() << " " << transf.translation().z());

        m_geoShapeConverter->decodeShape(clv->getShape());

        printChildren(cv);
    }
}

Trk::TrackingVolume* Muon::MuonInertMaterialBuilder::simplifyShape(Trk::TrackingVolume* trVol, bool blend) {
    // envelope
    const Trk::Volume* envelope = nullptr;
    // resolve composed volumes (returns constituents with material fraction accounting for subtractions & overlaps)
    std::vector<std::pair<const Trk::Volume*, std::pair<float, float> > > constituents = splitComposedVolume(trVol, m_simplify || blend);

    ATH_MSG_VERBOSE("simplifying shape for:" << trVol->volumeName() << "," << constituents[0].second.first << ","
                                             << constituents[0].second.second << " nr constituents " << constituents.size());
    for (unsigned int i = 0; i < constituents.size(); i++)
        ATH_MSG_VERBOSE("constituent:" << calculateVolume(constituents[i].first) << "," << constituents[i].second.first << ","
                                       << constituents[i].second.second);

    int simpleMode = 0;

    if (constituents.size() == 1) {  // easy case
        envelope = new Trk::Volume(*(constituents.front().first), trVol->transform());
        if (constituents.front().second.first > 0.999) simpleMode = 1;  // no need to simplify nor envelope
    } else {                                                            // construct envelope using constituent edges
        ATH_MSG_VERBOSE(" createEnvelope ");
        envelope = createEnvelope(trVol->transform(), constituents);
    }
    if (!envelope) envelope = new Trk::Volume(*trVol);

    // simplification

    Trk::TrackingVolume* newVol = nullptr;
    std::vector<const Trk::TrackingVolume*>* confinedVols = new std::vector<const Trk::TrackingVolume*>;

    std::string envName = trVol->volumeName();

    if (simpleMode == 1) {
        newVol = trVol;
        delete confinedVols;
    } else if (m_simplify && constituents.front().second.first > 0.5) {
        simpleMode = 2;
        if (constituents.size() == 1) {  // simplified volume : the scale factor refers to the density scaling
            double fraction = constituents.front().second.first;
            // simplified material rescales X0, l0 and density
            ATH_MSG_VERBOSE(" Applying scaling for " << trVol->volumeName() << " fraction " << fraction);

            Trk::Material mat(trVol->X0 / fraction, trVol->L0 / fraction, trVol->A, trVol->Z, fraction * trVol->rho);
            newVol = new Trk::TrackingVolume(*envelope, mat, nullptr, nullptr, envName);
            delete trVol;
            delete confinedVols;
        } else {  // enclose simplified constituents
            for (unsigned int ic = 0; ic < constituents.size(); ic++) {
                double fraction = constituents[ic].second.first;
                ATH_MSG_VERBOSE(" Applying scaling for " << trVol->volumeName() << " fraction " << fraction);
                // simplified material rescales X0, l0 and density
                Trk::Material mat(trVol->X0 / fraction, trVol->L0 / fraction, trVol->A, trVol->Z, fraction * trVol->rho);
                Trk::TrackingVolume* trc = new Trk::TrackingVolume(*(constituents[ic].first), mat, nullptr, nullptr, trVol->volumeName());
                confinedVols->push_back(trc);
            }
            envName = trVol->volumeName() + "_envelope";
            newVol = new Trk::TrackingVolume(*envelope, m_muonMaterial, confinedVols, envName);
            for (unsigned int iv = 0; iv < confinedVols->size(); iv++)
                Trk::TrackingVolumeManipulator::confineVolume(*((*confinedVols)[iv]), newVol);
            delete trVol;
        }
    } else {  // enclose the exact transcript
        ATH_MSG_VERBOSE(" enclose the exact transcript ");
        confinedVols->push_back(trVol);
        envName = trVol->volumeName() + "_envelope";
        newVol = new Trk::TrackingVolume(*envelope, m_muonMaterial, confinedVols, envName);
        Trk::TrackingVolumeManipulator::confineVolume(*trVol, newVol);
        //    for (unsigned int iv = 0; iv < confinedVols->size(); iv++)
        //	Trk::TrackingVolumeManipulator::confineVolume(*((*confinedVols)[iv]),newVol);
    }

    if (blend) {
        // save calculable volumes for blending
        std::vector<std::pair<std::unique_ptr<const Trk::Volume>, float> > confinedConst;
        for (unsigned int ic = 0; ic < constituents.size(); ic++) {
            float scale = simpleMode == 2 ? 1 : constituents[ic].second.first;
            confinedConst.emplace_back(
                std::make_unique<Trk::Volume>(*(constituents[ic].first), newVol->transform().inverse()), scale);
        }
        m_constituents.push_back(std::move(confinedConst));
    }

    delete envelope;
    return newVol;
}

double Muon::MuonInertMaterialBuilder::calculateVolume(const Trk::Volume* envelope) const {
    double envVol = 0.;

    if (!envelope) return 0.;

    const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*>(&(envelope->volumeBounds()));
    const Trk::CuboidVolumeBounds* box = dynamic_cast<const Trk::CuboidVolumeBounds*>(&(envelope->volumeBounds()));
    const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*>(&(envelope->volumeBounds()));
    const Trk::BevelledCylinderVolumeBounds* bcyl = dynamic_cast<const Trk::BevelledCylinderVolumeBounds*>(&(envelope->volumeBounds()));
    const Trk::PrismVolumeBounds* prism = dynamic_cast<const Trk::PrismVolumeBounds*>(&(envelope->volumeBounds()));
    const Trk::SimplePolygonBrepVolumeBounds* spb = dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*>(&(envelope->volumeBounds()));
    const Trk::CombinedVolumeBounds* comb = dynamic_cast<const Trk::CombinedVolumeBounds*>(&(envelope->volumeBounds()));
    const Trk::SubtractedVolumeBounds* sub = dynamic_cast<const Trk::SubtractedVolumeBounds*>(&(envelope->volumeBounds()));

    if (cyl)
        envVol = 2 * cyl->halfPhiSector() * (cyl->outerRadius() * cyl->outerRadius() - cyl->innerRadius() * cyl->innerRadius()) *
                 cyl->halflengthZ();
    if (box) envVol = (8 * box->halflengthX() * box->halflengthY() * box->halflengthZ());
    if (trd) envVol = (4 * (trd->minHalflengthX() + trd->maxHalflengthX()) * trd->halflengthY() * trd->halflengthZ());
    if (bcyl) {
        int type = bcyl->type();
        if (type < 1)
            envVol = 2 * bcyl->halfPhiSector() * (bcyl->outerRadius() * bcyl->outerRadius() - bcyl->innerRadius() * bcyl->innerRadius()) *
                     bcyl->halflengthZ();
        if (type == 1)
            envVol = 2 * bcyl->halflengthZ() *
                     (bcyl->halfPhiSector() * bcyl->outerRadius() * bcyl->outerRadius() -
                      bcyl->innerRadius() * bcyl->innerRadius() * tan(bcyl->halfPhiSector()));
        if (type == 2)
            envVol = 2 * bcyl->halflengthZ() *
                     (-bcyl->halfPhiSector() * bcyl->innerRadius() * bcyl->innerRadius() +
                      bcyl->outerRadius() * bcyl->outerRadius() * tan(bcyl->halfPhiSector()));
        if (type == 3)
            envVol = 2 * bcyl->halflengthZ() * tan(bcyl->halfPhiSector()) *
                     (bcyl->outerRadius() * bcyl->outerRadius() - bcyl->innerRadius() * bcyl->innerRadius());
    }
    if (prism) {
#ifdef TRKDETDESCR_USEFLOATPRECISON
#define double float
#endif
        std::vector<std::pair<double, double> > v = prism->xyVertices();
#ifdef TRKDETDESCR_USEFLOATPRECISON
#undef double
#endif
        double a2 = v[1].first * v[1].first + v[1].second * v[1].second + v[0].first * v[0].first + v[0].second * v[0].second -
                    2 * (v[0].first * v[1].first + v[0].second * v[1].second);
        double c2 = v[2].first * v[2].first + v[2].second * v[2].second + v[0].first * v[0].first + v[0].second * v[0].second -
                    2 * (v[0].first * v[2].first + v[0].second * v[2].second);
        double ca = v[1].first * v[2].first + v[1].second * v[2].second + v[0].first * v[0].first + v[0].second * v[0].second -
                    v[0].first * v[1].first - v[0].second * v[1].second - v[0].first * v[2].first - v[0].second * v[2].second;

        double vv2 = (a2 * c2 - ca * ca);
        double vv = 0.;
        // fix nans
        if (vv2 > 0.) vv = sqrt(vv2);
        envVol = vv * prism->halflengthZ();
        ATH_MSG_VERBOSE(" prism " << envVol << " vv " << vv << " a2*c2-ca*ca " << a2 * c2 - ca * ca << " v size " << v.size());
    }
    if (spb) {
        envVol = calculateVolume(spb->combinedVolume());  // exceptional use of combined volume (no intersections)
        ATH_MSG_VERBOSE(" spb volume " << envVol << " halflengthZ " << spb->halflengthZ());
    }
    if (comb) {
        ATH_MSG_VERBOSE(" comb volume first " << calculateVolume(comb->first()) << " second " << calculateVolume(comb->second()));
        envVol = calculateVolume(comb->first()) + calculateVolume(comb->second());
    }
    if (sub) { return -1; }

    return envVol;
}

void Muon::MuonInertMaterialBuilder::getObjsForTranslation(
    const GeoVPhysVol* pv, const Amg::Transform3D& transform,
    std::vector<std::pair<const GeoLogVol*, std::vector<Amg::Transform3D> > >& vols) const {
    // subcomponents
    unsigned int nc = pv->getNChildVols();
    ATH_MSG_VERBOSE(" INERT getObjsForTranslation from:" << pv->getLogVol()->getName() << "," << pv->getLogVol()->getMaterial()->getName()
                                                         << ", looping over " << nc << " children");
    for (unsigned int ic = 0; ic < nc; ic++) {
        Amg::Transform3D transf = pv->getXToChildVol(ic);
        const GeoVPhysVol* cv = &(*(pv->getChildVol(ic)));
        const GeoLogVol* clv = cv->getLogVol();
        if (!cv->getNChildVols()) {
            bool found = false;
            for (unsigned int is = 0; is < vols.size(); is++) {
                if (clv->getName() == vols[is].first->getName()) {
                    found = true;
                    vols[is].second.push_back(transform * transf);
                    break;
                }
            }
            if (!found) {
                std::vector<Amg::Transform3D> volTr;
                volTr.push_back(transform * transf);
                vols.emplace_back(clv, volTr);
                ATH_MSG_VERBOSE("INERT new volume added:" << clv->getName() << "," << clv->getMaterial()->getName());
                if (msg().level() <= MSG::VERBOSE) printInfo(cv);
            }
        } else {
            getObjsForTranslation(cv, transform * transf, vols);
        }
    }
}

const Trk::Volume* Muon::MuonInertMaterialBuilder::createEnvelope(
    const Amg::Transform3D& transf, std::vector<std::pair<const Trk::Volume*, std::pair<float, float> > > constituents) const {
    Trk::Volume* envelope = nullptr;

    std::vector<std::pair<const Trk::Volume*, std::pair<float, float> > >::iterator sIter = constituents.begin();
    std::vector<Amg::Vector3D> edges;
    bool cylinder = false;

    // double cVol = 0.;
    while (sIter != constituents.end()) {
        // cVol +=(*sIter).second.first * calculateVolume((*sIter).first);
        const Trk::SimplePolygonBrepVolumeBounds* spbBounds =
            dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*>(&((*sIter).first->volumeBounds()));
        const Trk::CylinderVolumeBounds* cylBounds = dynamic_cast<const Trk::CylinderVolumeBounds*>(&((*sIter).first->volumeBounds()));
        const Trk::CuboidVolumeBounds* cubBounds = dynamic_cast<const Trk::CuboidVolumeBounds*>(&((*sIter).first->volumeBounds()));
        const Trk::TrapezoidVolumeBounds* trdBounds = dynamic_cast<const Trk::TrapezoidVolumeBounds*>(&((*sIter).first->volumeBounds()));

        if (cylBounds) {
            double rOut = cylBounds->outerRadius();
            double hZ = cylBounds->halflengthZ();
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(rOut, rOut, hZ));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-rOut, rOut, hZ));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(rOut, -rOut, hZ));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-rOut, -rOut, hZ));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(rOut, rOut, -hZ));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-rOut, rOut, -hZ));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(rOut, -rOut, -hZ));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-rOut, -rOut, -hZ));
            cylinder = true;
            ATH_MSG_VERBOSE(" createEnvelope cylBounds ");
        } else if (cubBounds) {
            double x = cubBounds->halflengthX();
            double y = cubBounds->halflengthY();
            double z = cubBounds->halflengthZ();
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(x, y, z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-x, y, z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(x, -y, z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-x, -y, z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(x, y, -z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-x, y, -z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(x, -y, -z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-x, -y, -z));
            ATH_MSG_VERBOSE(" createEnvelope cubBounds ");
        } else if (trdBounds) {
            double x1 = trdBounds->minHalflengthX();
            double x2 = trdBounds->maxHalflengthX();
            double y = trdBounds->halflengthY();
            double z = trdBounds->halflengthZ();
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(x2, y, z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-x2, y, z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(x1, -y, z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-x1, -y, z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(x2, y, -z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-x2, y, -z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(x1, -y, -z));
            edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(-x1, -y, -z));
            ATH_MSG_VERBOSE(" createEnvelope trpBounds ");
        } else if (spbBounds) {
#ifdef TRKDETDESCR_USEFLOATPRECISON
#define double float
#endif
            std::vector<std::pair<double, double> > xyVtx = spbBounds->xyVertices();
#ifdef TRKDETDESCR_USEFLOATPRECISON
#undef double
#endif
            double z = spbBounds->halflengthZ();
            for (unsigned int iv = 0; iv < xyVtx.size(); iv++) {
                edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(xyVtx[iv].first, xyVtx[iv].second, -z));
                edges.push_back(transf.inverse() * (*sIter).first->transform() * Amg::Vector3D(xyVtx[iv].first, xyVtx[iv].second, z));
            }
            ATH_MSG_VERBOSE(" createEnvelope psbBounds ");
        } else {
            ATH_MSG_VERBOSE(" bounds not recognized ");
            return nullptr;
        }
        ++sIter;
    }

    double xMin = 13000.;
    double xMax = -13000.;
    double yMin = 13000.;
    double yMax = -13000.;
    double zMin = 25000.;
    double zMax = -25000.;
    // double rSize = -1.;
    for (unsigned int ie = 0; ie < edges.size(); ie++) {
        xMin = fmin(xMin, edges[ie].x());
        xMax = fmax(xMax, edges[ie].x());
        yMin = fmin(yMin, edges[ie].y());
        yMax = fmax(yMax, edges[ie].y());
        zMin = fmin(zMin, edges[ie].z());
        zMax = fmax(zMax, edges[ie].z());
        ATH_MSG_VERBOSE(" Envelope edges x " << edges[ie].x() << " y " << edges[ie].y() << " z " << edges[ie].z());
    }

    double xSize = 0.5 * (xMax - xMin);
    double ySize = 0.5 * (yMax - yMin);
    double zSize = 0.5 * (zMax - zMin);

    // bool cylEnv = false;

    if (cylinder && fabs(xSize - ySize) / fmax(xSize, ySize) < 0.1) {  // make it a cylinder
        envelope = new Trk::Volume(
            new Amg::Transform3D(transf * Amg::Translation3D(Amg::Vector3D(0.5 * (xMin + xMax), 0.5 * (yMin + yMax), 0.5 * (zMin + zMax)))),
            new Trk::CylinderVolumeBounds(sqrt(xSize * xSize + ySize * ySize), zSize));
        ATH_MSG_VERBOSE(" envelop cylinder radius "
                        << sqrt(xSize * xSize + ySize * ySize) << " zSize " << zSize << " xMin+xMax " << xMin + xMax << " yMin+yMax "
                        << yMin + yMax << " zMin+zMax " << zMin + zMax << " global z "
                        << ((transf * Amg::Translation3D(Amg::Vector3D(0.5 * (xMin + xMax), 0.5 * (yMin + yMax), 0.5 * (zMin + zMax)))) *
                            Amg::Vector3D(0., 0., 0.))
                               .z());

        // cylEnv = true;
    } else {
        envelope = new Trk::Volume(
            new Amg::Transform3D(transf * Amg::Translation3D(Amg::Vector3D(0.5 * (xMin + xMax), 0.5 * (yMin + yMax), 0.5 * (zMin + zMax)))),
            new Trk::CuboidVolumeBounds(xSize, ySize, zSize));

        ATH_MSG_VERBOSE(" envelop cube  xSize " << xSize << " ySize " << ySize << " zSize " << zSize << " xMin+xMax " << xMin + xMax
                                                << " yMin+yMax " << yMin + yMax << " zMin+zMax " << zMin + zMax);
    }

    /*
    // check if all edges really confined:
    bool envFine = true;
    for (unsigned int ie=0; ie < edges.size(); ie++) {
      if (!envelope->inside(edges[ie],0.001)) {
        std::cout << "something's wrong with the envelope, edge not confined"<< std::endl;
        envFine = false;
      }
    }
    if (!envFine) envelope = 0;
    */

    // std::cout << "volume/envelope fraction:"<< cVol/calculateVolume(envelope)<< std::endl;

    return envelope;
}

Amg::Vector3D Muon::MuonInertMaterialBuilder::getScanPoint(const Trk::Volume* vol) const {
    Amg::Vector3D gp(0., 0., 0.);

    const Trk::CuboidVolumeBounds* cub = dynamic_cast<const Trk::CuboidVolumeBounds*>(&(vol->volumeBounds()));
    const Trk::CylinderVolumeBounds* cyl = dynamic_cast<const Trk::CylinderVolumeBounds*>(&(vol->volumeBounds()));
    const Trk::TrapezoidVolumeBounds* trd = dynamic_cast<const Trk::TrapezoidVolumeBounds*>(&(vol->volumeBounds()));
    const Trk::SimplePolygonBrepVolumeBounds* spb = dynamic_cast<const Trk::SimplePolygonBrepVolumeBounds*>(&(vol->volumeBounds()));

    if (!cub && !cyl && !trd && !spb) {
        ATH_MSG_DEBUG("unknown scan shape!");
        return (vol->transform() * gp);
    }

    std::vector<double> rndm(3);
    // generate random numbers
    for (unsigned int ir = 0; ir < 3; ir++) rndm[ir] = m_flatDist->shoot();

    if (cub) {
        double x = cub->halflengthX();
        double y = cub->halflengthY();
        double z = cub->halflengthZ();
        gp = Amg::Vector3D(-x + 2 * x * rndm[0], -y + 2 * y * rndm[1], -z + 2 * z * rndm[2]);
    } else if (trd) {
        double x1 = trd->minHalflengthX();
        double x2 = trd->maxHalflengthX();
        double y = trd->halflengthY();
        double z = trd->halflengthZ();
        gp = Amg::Vector3D(-x2 + (x1 + x2) * rndm[0], -y + 2 * y * rndm[1], -z + 2 * z * rndm[2]);
        if (!vol->inside(vol->transform() * gp, 0.001)) {
            gp = Amg::Vector3D(x1 + (x1 + x2) * rndm[0], y - 2 * y * rndm[1], -z + 2 * z * rndm[2]);
            if (!vol->inside(vol->transform() * gp, 0.001)) ATH_MSG_DEBUG("trapezoid hit not correct:" << -y + 2 * y * rndm[1]);
        }
    } else if (cyl) {
        double x1 = cyl->innerRadius();
        double x2 = cyl->outerRadius();
        double y = cyl->halfPhiSector();
        double z = cyl->halflengthZ();
        double r = sqrt(x1 * x1 + rndm[0] * (x2 * x2 - x1 * x1));
        double phi = -y + 2 * y * rndm[1];
        gp = Amg::Vector3D(r * cos(phi), r * sin(phi), -z + 2 * z * rndm[2]);
    } else if (spb) {
        double z = spb->halflengthZ();
        std::vector<const Trk::Volume*> subVols;
        std::vector<double> subVolumes;
        const Trk::Volume* comVol = spb->combinedVolume();
        const Trk::CombinedVolumeBounds* comb = dynamic_cast<const Trk::CombinedVolumeBounds*>(&(comVol->volumeBounds()));
        while (comb) {
            const Trk::CombinedVolumeBounds* comb1 = dynamic_cast<const Trk::CombinedVolumeBounds*>(&(comb->first()->volumeBounds()));
            const Trk::CombinedVolumeBounds* comb2 = dynamic_cast<const Trk::CombinedVolumeBounds*>(&(comb->second()->volumeBounds()));
            if (!comb1 && !comb2) {
                subVols.push_back(comb->second());
                comVol = comb->first();
                comb = nullptr;
            } else if (comb1) {
                comb = comb1;
                subVols.push_back(comb->second());
            } else {
                comb = comb2;
                subVols.push_back(comb->first());
            }
        }
        subVols.push_back(comVol);
        subVolumes.resize(subVols.size());
        double spbVol = 0.;
        for (unsigned int i = 0; i < subVols.size(); i++) {
            subVolumes[i] = calculateVolume(subVols[i]);
            spbVol += subVolumes[i];
        }
        // first rndm defines prism
        unsigned int iPrism = 0;
        double vCount = 0.;
        while (iPrism < subVols.size() && vCount < rndm[0] * spbVol) {
            vCount += subVolumes[iPrism];
            iPrism++;
        }
        vCount = subVolumes[iPrism - 1] - (vCount - rndm[0] * spbVol);
        double zfr = vCount / subVolumes[iPrism - 1];
        // std::cout << "prism picked:"<< iPrism-1<<","<< subVolumes.size()<<": z fraction:"<<zfr<< std::endl;
        // xy
        const Trk::PrismVolumeBounds* prism = dynamic_cast<const Trk::PrismVolumeBounds*>(&(subVols[iPrism - 1]->volumeBounds()));
        if (prism) {
#ifdef TRKDETDESCR_USEFLOATPRECISON
#define double float
#endif
            std::vector<std::pair<double, double> > xy = prism->xyVertices();
#ifdef TRKDETDESCR_USEFLOATPRECISON
#undef double
#endif
            // std::cout << "prism vertices:"<< xy.size()<< std::endl;
            gp = Amg::Vector3D(xy[2].first + sqrt(rndm[1]) * (xy[0].first - xy[2].first + rndm[2] * (xy[1].first - xy[0].first)),
                               xy[2].second + sqrt(rndm[1]) * (xy[0].second - xy[2].second + rndm[2] * (xy[1].second - xy[0].second)),
                               -z + zfr * 2 * z);
        }
    }

    if (!vol->inside(vol->transform() * gp, 0.001)) ATH_MSG_DEBUG("test hit:wrong scan hit:" << gp);

    return (vol->transform() * gp);
}

std::vector<std::pair<const Trk::Volume*, std::pair<float, float> > > Muon::MuonInertMaterialBuilder::splitComposedVolume(
    const Trk::Volume* trVol, bool estimateVol) const {
    std::vector<const Trk::Volume*> garbage;
    std::vector<std::pair<const Trk::Volume*, const Trk::Volume*> > constituents;
    constituents.emplace_back(trVol, nullptr);
    std::vector<std::pair<const Trk::Volume*, const Trk::Volume*> >::iterator sIter = constituents.begin();
    const Trk::Volume* subVol = nullptr;
    while (sIter != constituents.end()) {
        const Trk::CombinedVolumeBounds* comb = dynamic_cast<const Trk::CombinedVolumeBounds*>(&((*sIter).first->volumeBounds()));
        const Trk::SubtractedVolumeBounds* sub = dynamic_cast<const Trk::SubtractedVolumeBounds*>(&((*sIter).first->volumeBounds()));
        if (comb) {
            subVol = (*sIter).second;
            sIter = constituents.erase(sIter);
            if (comb->intersection()) {
                Trk::Volume* newSubVol =
                    new Trk::Volume(nullptr, new Trk::SubtractedVolumeBounds(comb->first()->clone(), comb->second()->clone()));
                if (subVol) {
                    Trk::Volume* newCSubVol = new Trk::Volume(nullptr, new Trk::CombinedVolumeBounds(subVol->clone(), newSubVol, false));
                    constituents.insert(sIter, std::pair<const Trk::Volume*, const Trk::Volume*>(comb->first(), newCSubVol));
                    garbage.push_back(newCSubVol);
                } else {
                    constituents.insert(sIter, std::pair<const Trk::Volume*, const Trk::Volume*>(comb->first(), newSubVol));
                    garbage.push_back(newSubVol);
                }
            } else {
                constituents.insert(sIter, std::pair<const Trk::Volume*, const Trk::Volume*>(comb->first(), subVol));
                constituents.insert(sIter, std::pair<const Trk::Volume*, const Trk::Volume*>(comb->second(), subVol));
            }
            sIter = constituents.begin();
        } else if (sub) {
            subVol = (*sIter).second;
            sIter = constituents.erase(sIter);
            if (subVol) {
                Trk::Volume* newSubVol = new Trk::Volume(nullptr, new Trk::CombinedVolumeBounds(subVol->clone(), sub->inner()->clone(), false));
                constituents.insert(sIter, std::pair<const Trk::Volume*, const Trk::Volume*>(sub->outer(), newSubVol));
                garbage.push_back(newSubVol);
            } else {
                constituents.insert(sIter, std::pair<const Trk::Volume*, const Trk::Volume*>(sub->outer(), sub->inner()));
            }
            sIter = constituents.begin();
        } else {
            ++sIter;
        }
    }

    // estimate fraction of cutouts & overlaps
    std::vector<std::pair<const Trk::Volume*, std::pair<float, float> > > wConst;
    for (unsigned int i = 0; i < constituents.size(); i++)
        wConst.emplace_back(constituents[i].first, std::pair<float, float>(1., 0.));

    if (estimateVol && (constituents.size() > 1 || constituents[0].second)) {
        for (unsigned int iv = 0; iv < constituents.size(); iv++) {
            const Trk::Volume* replaceVol = nullptr;
            // check if volume can be optimized (replace cylinder by tube)
            if (constituents[iv].second) {
                const Trk::CylinderVolumeBounds* cyl =
                    dynamic_cast<const Trk::CylinderVolumeBounds*>(&(constituents[iv].first->volumeBounds()));
                if (cyl && cyl->innerRadius() < 0.001) {
                    // recalculate inner radius for better efficiency
                    std::vector<std::pair<const Trk::Volume*, std::pair<float, float> > > subtr =
                        splitComposedVolume(constituents[iv].second, false);
                    for (unsigned int is = 0; is < subtr.size(); is++) {
                        const Trk::CylinderVolumeBounds* cyls =
                            dynamic_cast<const Trk::CylinderVolumeBounds*>(&(subtr[is].first->volumeBounds()));
                        if (cyls) {
                            Amg::Vector3D dc = subtr[is].first->transform().inverse() * constituents[iv].first->center();
                            if (cyls->outerRadius() < cyl->outerRadius() && dc.perp() < 0.001 && fabs(dc.z()) < 10. &&
                                fabs(cyl->halflengthZ() - cyls->halflengthZ()) < 10.) {
                                replaceVol = new Trk::Volume(new Amg::Transform3D(constituents[iv].first->transform()),
                                                             new Trk::CylinderVolumeBounds(cyls->outerRadius(), cyl->outerRadius(),
                                                                                           cyl->halfPhiSector(), cyl->halflengthZ()));
                                break;
                            }
                        }
                    }
                }
            }
            // hit & miss sampling
            float cutout = 0.;
            float overlap = 0.;
            const Trk::Volume* scanVol = replaceVol ? replaceVol : constituents[iv].first;
            double nHits = calculateVolume(scanVol) / 1.e6;
            if (nHits < 100.) nHits = 100.;
            if (nHits > 1000.) nHits = 1000.;
            int cut = 1;
            int ovl = 0;
            for (unsigned int ih = 0; ih < nHits; ih++) {
                Amg::Vector3D gp = getScanPoint(scanVol);
                cut = 1;
                ovl = 0;
                if (constituents[iv].second && constituents[iv].second->inside(gp, 0.))
                    cut = 0;
                else {
                    for (unsigned int icv = 0; icv < constituents.size(); icv++)
                        if (icv != iv && constituents[icv].first->inside(gp, 0.) &&
                            (!constituents[icv].second || !constituents[icv].second->inside(gp, 0.)))
                            ovl = 1;
                }
                cutout += cut;
                overlap += ovl;
            }
            if (replaceVol) wConst[iv].first = replaceVol;
            wConst[iv].second.first = cutout / nHits;
            wConst[iv].second.second = overlap / nHits;
            // std::cout << " cutout " << cutout/nHits << " overlap " << overlap/nHits << std::endl;
        }
    }

    for (size_t i = 0; i < garbage.size(); i++) delete garbage[i];

    return wConst;
}
