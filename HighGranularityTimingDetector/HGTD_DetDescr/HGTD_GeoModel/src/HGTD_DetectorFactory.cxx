/*
 Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 */

/*  Three-ring detector layout, created by Christian and David
 *  Sep 2020 Lianyou SHAN merged
 *  hgtdlayout-master/place_modules_option2.py by Christina Agapopoulou
 *  into LArGeoModel/LArGeoEndcap/EndcapCryostatConstruction in R20.20.14.6
 *  Dec 2020 Lianyou SHAN integrate 2-ring layout
 */

#include "HGTD_GeoModel/HGTD_DetectorFactory.h"

#include <string>
#include <sstream>
#include <math.h>

#include "AthenaBaseComps/AthMsgStreamMacros.h"
#include "AthenaBaseComps/AthCheckMacros.h"
#include "StoreGate/StoreGateSvc.h" // For alignment getAlignableFolderType()
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "DetDescrConditions/AlignableTransformContainer.h"

#include "HGTD_GeoModel/HGTD_GeoModelAthenaComps.h"
#include "GeoModelKernel/GeoNameTag.h"
#include "GeoModelKernel/GeoIdentifierTag.h"
#include "GeoModelKernel/GeoMaterial.h"
#include "GeoModelKernel/GeoTransform.h"
#include "GeoModelKernel/GeoAlignableTransform.h"
#include "GeoModelKernel/GeoTorus.h"
#include "GeoModelKernel/GeoTube.h"
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "GeoModelInterfaces/IGeoModelSvc.h"
#include "GeoModelInterfaces/StoredMaterialManager.h"
#include "GeoModelUtilities/DecodeVersionKey.h"
#include "GeoModelUtilities/StoredPhysVol.h"
#include "GeoModelUtilities/StoredAlignX.h"
#include "GaudiKernel/SystemOfUnits.h"

#include "HGTD_Identifier/HGTD_ID.h"

#include "HGTD_ReadoutGeometry/HGTD_DetectorManager.h"
#include "HGTD_ReadoutGeometry/HGTD_ModuleDesign.h"
#include "HGTD_ReadoutGeometry/HGTD_DetectorElement.h"

#include "ReadoutGeometryBase/PixelDiodeMatrix.h"
#include "ReadoutGeometryBase/SiCommonItems.h"

#include "RDBAccessSvc/IRDBAccessSvc.h"
#include "RDBAccessSvc/IRDBRecordset.h"
#include "RDBAccessSvc/IRDBRecord.h"

using namespace std;
using namespace InDetDD;

HGTD_DetectorFactory::HGTD_DetectorFactory( HGTD_GeoModelAthenaComps* athComps ) :
  InDetDD::DetectorFactoryBase( athComps ),
  m_athComps( athComps ),
  m_materialMgr( nullptr ),
  m_geomVersion( -1 ),
  m_outputIdfr( false ) {
    // create the detector manager
    m_detectorManager = new HGTD_DetectorManager( detStore() );

    ATH_MSG_INFO( "HGTD geometry from hard-coded definition - No Information being taken from Geometry Tag!" );

    
    // Create SiCommonItems. These are items that are shared by all elements
    m_commonItems = std::make_unique<const InDetDD::SiCommonItems>(m_athComps->getIdHelper());
    
    // temporarily hardcode the HGTD version to build until the geo db has been updated with tables for 3-ring layout
    // m_geomVersion = 0; // two-ring layout
    m_geomVersion = 1; // three-ring layout
}

HGTD_DetectorFactory::~HGTD_DetectorFactory() {
    // NB the detector manager (m_detectorManager) is stored in the detector store by the Tool and so we don't delete it.
}

void HGTD_DetectorFactory::setPrintIdentifierDict( bool print ) {
    m_outputIdfr = print;
    return;
}

void HGTD_DetectorFactory::create(GeoPhysVol* world) {

    ATH_MSG_INFO( "Building HGTD detector");

    // initialize the geometry parameters
    // Originally taken from geometry db, but now fully hard-coded until move to GMX implementation
    initializeGeoParameters();

    // for now the position of the HGTD mother volumes is hardcoded - TODO: take from db!
    float zMother = 3482.5;

    // build logical volumes for the two endcaps
    const GeoLogVol* positiveEndcapLogicalVolume = buildEndcapLogicalVolume(true);
    const GeoLogVol* negativeEndcapLogicalVolume = buildEndcapLogicalVolume(false);

    // create and place positive endcap
    world->add(new GeoNameTag("HGTD_Pos"));
    world->add(new GeoIdentifierTag(9));
    ATH_MSG_INFO( "HGTD_Pos mother volume will be placed at z = " << zMother << " mm" );
    world->add(new GeoTransform(GeoTrf::TranslateZ3D(zMother)));
    GeoVPhysVol* endcapPos = build( positiveEndcapLogicalVolume,  true);
    world->add( endcapPos );
    m_detectorManager->addTreeTop( endcapPos);

    // create and place negative endcap
    world->add(new GeoNameTag("HGTD_Neg"));
    world->add(new GeoIdentifierTag(-9));
    ATH_MSG_INFO( "HGTD_Neg mother volume will be placed at z = " << -zMother << " mm" );
    world->add(new GeoTransform(GeoTrf::TranslateZ3D(-zMother)));
    world->add(new GeoTransform(GeoTrf::RotateY3D(180.0*Gaudi::Units::deg)));
    GeoVPhysVol* endcapNeg = build( negativeEndcapLogicalVolume, false);
    world->add( endcapNeg );
    m_detectorManager->addTreeTop( endcapNeg );
    
    // Add SiCommonItems to HGTD_DetectorManager to hold and delete it.
    m_detectorManager->setCommonItems(std::move(m_commonItems));

    return;
}

// initialize the geometry parameters (now fully hard-coded until move to GMX implementation)
void HGTD_DetectorFactory::initializeGeoParameters() {

    // retrieve the material manager (can't use ATH_CHECK macros within create(), it seems..)
    StatusCode sc = detStore()->retrieve(m_materialMgr, std::string("MATERIALS"));
    if (sc != StatusCode::SUCCESS) {
      ATH_MSG_ERROR("Cannot retrieve material manager from DetStore");
    }

    // temporarily hard-code custom materials - eventually will be defined in xml once HGTD migrates to GeoModelXML detector description

    GeoMaterial* CFRP = new GeoMaterial("hgtd::CFRP", 1.78*(CLHEP::gram / CLHEP::cm3)); // copy of sct::CFRP used in 21.9
    CFRP->add(m_materialMgr->getElement("Carbon"), 0.92);
    CFRP->add(m_materialMgr->getElement("Hydrogen"), 0.02);
    CFRP->add(m_materialMgr->getElement("Oxygen"), 0.05);
    m_materialMgr->addMaterial("hgtd", CFRP);

    GeoMaterial* Honeycomb = new GeoMaterial("hgtd::Honeycomb", 0.042*(CLHEP::gram / CLHEP::cm3)); // copy of muo::Honeycomb used in 21.9
    Honeycomb->add(m_materialMgr->getElement("Carbon"), 0.88);
    Honeycomb->add(m_materialMgr->getElement("Hydrogen"), 0.11);
    m_materialMgr->addMaterial("hgtd", Honeycomb);

    GeoMaterial* Peek = new GeoMaterial("hgtd::Peek", 1.3*(CLHEP::gram / CLHEP::cm3)); // copy of pix::Peek used in 21.9
    Peek->add(m_materialMgr->getElement("Hydrogen"), 0.04);
    Peek->add(m_materialMgr->getElement("Carbon"), 0.79);
    Peek->add(m_materialMgr->getElement("Oxygen"), 0.16);
    m_materialMgr->addMaterial("hgtd", Peek);

    GeoMaterial* CO2_Liquid = new GeoMaterial("hgtd::CO2_Liquid", 1.032*(CLHEP::gram / CLHEP::cm3)); // copy of pix::CO2_Liquid used in 21.9
    CO2_Liquid->add(m_materialMgr->getElement("Carbon"), 0.27);
    CO2_Liquid->add(m_materialMgr->getElement("Oxygen"), 0.72);
    m_materialMgr->addMaterial("hgtd", CO2_Liquid);

    GeoMaterial* CO2 = new GeoMaterial("hgtd::CO2", 0.001842*(CLHEP::gram / CLHEP::cm3)); // copy of trt::CO2 used in 21.9
    CO2->add(m_materialMgr->getElement("Carbon"), 0.27);
    CO2->add(m_materialMgr->getElement("Oxygen"), 0.72);
    m_materialMgr->addMaterial("hgtd", CO2);

    GeoMaterial* CFiberSupport = new GeoMaterial("hgtd::CFiberSupport", 0.189*(CLHEP::gram / CLHEP::cm3)); // copy of sct::CFiberSupport used in 21.9
    CFiberSupport->add(m_materialMgr->getElement("Carbon"), 1.0);
    m_materialMgr->addMaterial("hgtd", CFiberSupport);

    GeoMaterial* CuKapton = new GeoMaterial("hgtd::CuKapton", 2.94*(CLHEP::gram / CLHEP::cm3)); // copy of sct::CuKapton used in 21.9
    CuKapton->add(m_materialMgr->getElement("Copper"), 0.61);
    CuKapton->add(m_materialMgr->getElement("Carbon"), 0.26);
    CuKapton->add(m_materialMgr->getElement("Hydrogen"), 0.01);
    CuKapton->add(m_materialMgr->getElement("Oxygen"), 0.08);
    CuKapton->add(m_materialMgr->getElement("Nitrogen"), 0.02);
    m_materialMgr->addMaterial("hgtd", CuKapton);

    GeoMaterial* BoratedPolyethelyne = new GeoMaterial("hgtd::BoratedPolyethelyne", 0.99*(CLHEP::gram / CLHEP::cm3)); // copy of LAr::BoratedPolyethelyne used in 21.9
    BoratedPolyethelyne->add(m_materialMgr->getElement("Hydrogen"), 0.13);
    BoratedPolyethelyne->add(m_materialMgr->getElement("Carbon"), 0.81);
    BoratedPolyethelyne->add(m_materialMgr->getElement("Boron"), 0.05);
    m_materialMgr->addMaterial("hgtd", BoratedPolyethelyne);

    GeoMaterial* FEBoards = new GeoMaterial("hgtd::FEBoards", 0.99*(CLHEP::gram / CLHEP::cm3)); // copy of LAr::FEBoards used in 21.9
    FEBoards->add(m_materialMgr->getElement("Silicon"), 0.27);
    FEBoards->add(m_materialMgr->getElement("Oxygen"), 0.34);
    FEBoards->add(m_materialMgr->getElement("Copper"), 0.28);
    FEBoards->add(m_materialMgr->getElement("Hydrogen"), 0.01);
    FEBoards->add(m_materialMgr->getElement("Carbon"), 0.09);
    m_materialMgr->addMaterial("hgtd", FEBoards);

    GeoMaterial* Epoxy = new GeoMaterial("hgtd::Epoxy", 1*(CLHEP::gram / CLHEP::cm3)); // copy of sct::Epoxy used in 21.9
    Epoxy->add(m_materialMgr->getElement("Carbon"), 0.76);
    Epoxy->add(m_materialMgr->getElement("Hydrogen"), 0.07);
    Epoxy->add(m_materialMgr->getElement("Oxygen"), 0.16);
    m_materialMgr->addMaterial("hgtd", Epoxy);

    // Hardcoded box parameters taken from geometry db:

    // Node: HGTDBox                      BOX                   DX       DY       DZ       ZPOS     MATERIAL
    m_boxVolPars["HGTDModule0"]        = {"HGTDModule0",        11,      20,      1.75,    0,       "std::Air"};
    m_boxVolPars["HGTDModule1"]        = {"HGTDModule1",        11,      20,      1.75,    0,       "std::Air"};
    m_boxVolPars["HGTDModule2"]        = {"HGTDModule2",        11,      20,      1.75,    0,       "std::Air"};
    m_boxVolPars["HGTDModule3"]        = {"HGTDModule3",        11,      20,      1.75,    0,       "std::Air"};
    m_boxVolPars["HGTD::Hybrid"]       = {"HGTD::Hybrid",       10.25,   20,      .175,    0,       "hgtd::CuKapton"};
    m_boxVolPars["HGTD::GlueSensor"]   = {"HGTD::GlueSensor",   10.25,   20,      .04,     0,       "hgtd::Epoxy"};
    m_boxVolPars["HGTD::GlueAsic"]     = {"HGTD::GlueAsic",     11,      20,      .04,     0,       "hgtd::Epoxy"};
    m_boxVolPars["HGTDSiSensor0"]      = {"HGTDSiSensor0",      10.25,   20,      .025,    0,       "std::Silicon"};
    m_boxVolPars["HGTDSiSensor1"]      = {"HGTDSiSensor1",      10.25,   20,      .025,    0,       "std::Silicon"};
    m_boxVolPars["HGTDSiSensor2"]      = {"HGTDSiSensor2",      10.25,   20,      .025,    0,       "std::Silicon"};
    m_boxVolPars["HGTDSiSensor3"]      = {"HGTDSiSensor3",      10.25,   20,      .025,    0,       "std::Silicon"};
    m_boxVolPars["HGTD::LGADInactive"] = {"HGTD::LGADInactive", 10.25,   20,      .1,      0,       "std::Silicon"};
    m_boxVolPars["HGTD::ASIC"]         = {"HGTD::ASIC",         11,      20,      .15,     0,       "std::Silicon"};

    // Add a dummy entry that will be used to leave some space - no volume will actually be created for this
    // needed after fix of ASIC thickness (and material) in HGTD-TDR-01 tag (ATLAS-P2-ITK-17-04-02 and later), compared to HGTD-TDR-00
    double moduleSpaceHalfZ = 0.225;
    m_boxVolPars["HGTD::ModuleSpace"] = {"HGTD::ModuleSpace", 11, 20, moduleSpaceHalfZ, 0, "std::Air"};



    // Node: HGTDTubs                       TUBE                       RMIN     RMAX     DZ       ZPOS     MATERIAL
    m_cylVolPars["HGTD_mother"]          = {"HGTD_mother",             100,     1100,    62.5,    -3252,   "std::Air"};
    m_cylVolPars["HGTD::FrontCover"]     = {"HGTD::FrontCover",        120,     1000,    7.5,     0,       "hgtd::CFiberSupport"};
    m_cylVolPars["HGTD::FlexPackage"]    = {"HGTD::FlexPackage",       120,     660,     2,       0,       "std::Air"};
    m_cylVolPars["HGTD::FlexTube"]       = {"HGTD::FlexTube",          120,     660,     .175,    0,       "hgtd::CuKapton"};
    m_cylVolPars["HGTD::ModuleLayer0"]   = {"HGTD::ModuleLayer0",      120,     660,     3.75,    0,       "std::Air"};
    m_cylVolPars["HGTD::ModuleLayer1"]   = {"HGTD::ModuleLayer1",      120,     660,     3.75,    0,       "std::Air"};
    m_cylVolPars["HGTD::ModuleLayer2"]   = {"HGTD::ModuleLayer2",      120,     660,     3.75,    0,       "std::Air"};
    m_cylVolPars["HGTD::ModuleLayer3"]   = {"HGTD::ModuleLayer3",      120,     660,     3.75,    0,       "std::Air"};
    m_cylVolPars["HGTD::CoolingPlate"]   = {"HGTD::CoolingPlate",      120,     920,     3,       0,       "std::Aluminium"};
    m_cylVolPars["HGTD::SupportPlate"]   = {"HGTD::SupportPlate",      120,     660,     .5,      0,       "std::Aluminium"};
    m_cylVolPars["HGTD::ToleranceFront"] = {"HGTD::ToleranceFront",    120,     660,     1,       0,       "std::Air"};
    m_cylVolPars["HGTD::ToleranceBack"]  = {"HGTD::ToleranceBack",     120,     660,     1,       0,       "std::Air"};
    m_cylVolPars["HGTD::ToleranceMid"]   = {"HGTD::ToleranceMid",      120,     660,     1,       0,       "std::Air"};
    m_cylVolPars["HGTD::ModeratorIn"]    = {"HGTD::ModeratorIn",       120,     900,     15,      0,       "hgtd::BoratedPolyethelyne"};
    m_cylVolPars["HGTD::ModeratorOut"]   = {"HGTD::ModeratorOut",      120,     1100,    10,      0,       "hgtd::BoratedPolyethelyne"};
    m_cylVolPars["HGTD::BackCover"]      = {"HGTD::BackCover",         120,     1100,    4,       0,       "hgtd::CFiberSupport"};
    m_cylVolPars["HGTD::PeriphElec"]     = {"HGTD::PeriphElec",        674,     900,     1,       2,       "hgtd::FEBoards"};

    m_cylVolPars["HGTD::InnerRCover1"] = {"HGTD::InnerRCover1", 110., 111., 105./2, -10., "hgtd::CFRP"};
    // the InnerRCover bulk should be 70% aerogel and 30% honeycomb made from "aradime" (not defined - using "hgtd::Honeycomb" for now)
    // proportions should be 50/50 by weight, which is used for GeoMaterial fractions
    // TODO: these should be double-checked, or at least that the density/weight matches engineering drawings
    GeoMaterial* innerRCoverBulkMaterial = new GeoMaterial("hgtd::AerogelAndHoneycomb", 0.17*(CLHEP::gram / CLHEP::cm3));
    innerRCoverBulkMaterial->add(m_materialMgr->getMaterial("std::Aerogel"), 0.5);
    innerRCoverBulkMaterial->add(m_materialMgr->getMaterial("hgtd::Honeycomb"), 0.5);
    m_materialMgr->addMaterial("hgtd", innerRCoverBulkMaterial);
    m_cylVolPars["HGTD::InnerRCover2"] = {"HGTD::InnerRCover2", 111., 119., 105./2, -10., "hgtd::AerogelAndHoneycomb"};
    m_cylVolPars["HGTD::InnerRCover3"] = {"HGTD::InnerRCover3", 119., 120., 105./2, -10., "hgtd::CFRP"};
    m_cylVolPars["HGTD::OuterRCover"]  = {"HGTD::OuterRCover", 980., 1000., 82./2, -6.5, "hgtd::Peek"};
    m_cylVolPars["HGTD::PeripheralCoolingLines"] = {"HGTD::PeripheralCoolingLines", 920., 980., 3./2, 31., "std::SSteel"};
    // TODO: outer cover should be 40% "hgtd::Peek" and 60% electrical connectors (unclear material)

    m_cylVolPars["HGTD::CoolingTube"] = {"HGTD::CoolingTubes", 0, 0, 2.0, 0, "std::Titanium"};
    // Coolant should be 50% liquid and 50% gas CO2 ("hgtd::CO2")
    GeoMaterial* coolantMaterial = new GeoMaterial("hgtd::CO2CoolantMix", 0.55*(CLHEP::gram / CLHEP::cm3));
    coolantMaterial->add(m_materialMgr->getMaterial("hgtd::CO2_Liquid"), 0.5);
    coolantMaterial->add(m_materialMgr->getMaterial("hgtd::CO2"), 0.5);
    m_materialMgr->addMaterial("hgtd", coolantMaterial);
    m_cylVolPars["HGTD::CoolingTubeFluid"] = {"HGTD::CoolingTubeFluid", 0, 0, 1.5, 0, "hgtd::CO2CoolantMix"};



    // These parameters were not in the db (they don't fit into the cylinder or box structures used above)
    m_hgtdPars = { 320., // rMid
        640., // rOuter - only used in one place, and there 20 mm is added to it...
        0.,   // disk1Rotation (in degrees)
        15., // disk2Rotation (in degrees)
        1.,   // rowSpaceSide
        4.,   // rowBacksideInnerShift
        17.,  // rowBacksideOuterShift
        1.5,  // moduleSpaceInner
        12.5, // moduleSpaceOuter
        0.456 // flexSheetSpacing
    };

    return;
}

// prepare an envelope for one HGTD side
GeoLogVol* HGTD_DetectorFactory::buildEndcapLogicalVolume(bool isPositiveSide) {

    // build the solid volume
    GeoTube* world_solid_hgtd = new GeoTube(m_cylVolPars["HGTD_mother"].rMin, m_cylVolPars["HGTD_mother"].rMax,
                                            m_cylVolPars["HGTD_mother"].zHalf);

    // build the logical volume
    std::string name = isPositiveSide ? "HGTD_PositiveEndcap" : "HGTD_NegativeEndcap";
    GeoLogVol* world_logical_hgtd  = new GeoLogVol( name.c_str(), world_solid_hgtd,
                            m_materialMgr->getMaterial( m_cylVolPars[ "HGTD_mother"].material) );

    return world_logical_hgtd;
}


GeoVPhysVol* HGTD_DetectorFactory::build( const GeoLogVol* logicalEnvelope, bool bPos) {

    ATH_MSG_INFO( "**************************************************");
    ATH_MSG_INFO( "       Building HGTD geometry , side =  " << bPos << "    ");
    ATH_MSG_INFO( "**************************************************" );

    GeoFullPhysVol* HGTDparent = new GeoFullPhysVol( logicalEnvelope );

    // to be calculated from parameters in db using map
    double motherHalfZ = ((GeoTube*) HGTDparent->getLogVol()->getShape())->getZHalfLength();
    double modulePackageHalfZtot = 3.5/2 + 4./2; // including flex - can we not get this from the db numbers? /CO

    double modulePackageHalfZ = 2*m_boxVolPars["HGTD::GlueSensor"].zHalf + m_boxVolPars["HGTDSiSensor0"].zHalf
    + m_boxVolPars["HGTD::LGADInactive"].zHalf + m_boxVolPars["HGTD::ASIC"].zHalf
    + m_boxVolPars["HGTD::Hybrid"].zHalf + m_boxVolPars["HGTD::ModuleSpace"].zHalf;

    // add volumes by key name to ordered vector, outside in (from larger z to smaller)
    std::vector<std::string> hgtdVolumes;
    hgtdVolumes.push_back("HGTD::ModeratorOut"); // Out as in outside the vessel
    hgtdVolumes.push_back("HGTD::BackCover");
    hgtdVolumes.push_back("HGTD::ToleranceBack");
    hgtdVolumes.push_back("HGTD::ModeratorIn"); // In as in inside the vessel

    hgtdVolumes.push_back("HGTD::ModuleLayer3");
    hgtdVolumes.push_back("HGTD::SupportPlate");
    hgtdVolumes.push_back("HGTD::CoolingPlate");
    hgtdVolumes.push_back("HGTD::SupportPlate");
    hgtdVolumes.push_back("HGTD::ModuleLayer2");

    hgtdVolumes.push_back("HGTD::ToleranceMid");

    hgtdVolumes.push_back("HGTD::ModuleLayer1");
    hgtdVolumes.push_back("HGTD::SupportPlate");
    hgtdVolumes.push_back("HGTD::CoolingPlate");
    hgtdVolumes.push_back("HGTD::SupportPlate");
    hgtdVolumes.push_back("HGTD::ModuleLayer0");

    hgtdVolumes.push_back("HGTD::ToleranceFront");
    hgtdVolumes.push_back("HGTD::FrontCover");
    // Important - these must come last since they will otherwise shift positions of the previous volumes!
    hgtdVolumes.push_back("HGTD::InnerRCover1"); // don't reorder!
    hgtdVolumes.push_back("HGTD::InnerRCover2"); // don't reorder!
    hgtdVolumes.push_back("HGTD::InnerRCover3"); // don't reorder!
    hgtdVolumes.push_back("HGTD::OuterRCover"); // don't reorder!
    hgtdVolumes.push_back("HGTD::PeripheralCoolingLines"); // don't reorder!

    // Now build up the solid, logical and physical volumes as appropriate (starting from the outermost volume)
    // We first start with the volumes we'll reuse several times

    //////////////////////////
    // FLEX PACKAGE VOLUMES //
    //////////////////////////

    // Flex package volume modeled as 8 concentric flex sheets with progressively larger inner radius
    // Order of sheets depend on whether package is for front or back of a cooling plate
    // First calculate the inner radii for the flex sheets
    GeoCylVolParams packagePars = m_cylVolPars["HGTD::FlexPackage"];
    GeoCylVolParams flexPars = m_cylVolPars["HGTD::FlexTube"];
    std::vector<double> flexSheetInnerR;
    double currentInnerR = 144.; // adding flex sheets from the second sensor (all have the hybrid already)
    for (int flexSheet = 0; flexSheet < 8; flexSheet++) {
        flexSheetInnerR.push_back(currentInnerR);
        // set the inner radius for the next flex sheet, increased by two module heights and two radius-dependent spaces per sheet
        currentInnerR += m_boxVolPars["HGTDModule0"].xHalf*2 * (2 + 2 * (flexSheet < 4 ? 0.2 : 0.8) );
    }

    // build up the two flex volumes for front (0) and back (1) sides
    GeoPhysVol* flexPackagePhysical[2] = {};
    for (int flexVolume = 0; flexVolume < 2; flexVolume++) {
        std::vector<double> rInner = flexSheetInnerR;
        if (flexVolume) reverse(rInner.begin(), rInner.end()); // reverse order for backside flex package

        GeoTube*    flexPackageSolid = new GeoTube(packagePars.rMin, packagePars.rMax, packagePars.zHalf);
        GeoLogVol*  flexPackageLogical = new GeoLogVol(packagePars.name, flexPackageSolid, m_materialMgr->getMaterial(packagePars.material));
        flexPackagePhysical[flexVolume] = new GeoPhysVol(flexPackageLogical);
        // build up a volume of flex cables, starting in z at half a flex layer from the edge of the flex package volume
        double flexZoffset = packagePars.zHalf - flexPars.zHalf;
        for (int flexSheet = 0; flexSheet < 8; flexSheet++) {
            GeoTube*    hgtdFlexSolid    = new GeoTube(rInner[flexSheet], flexPars.rMax, flexPars.zHalf);
            GeoLogVol*  hgtdFlexLogical  = new GeoLogVol("HGTD::FlexTube"+std::to_string(flexSheet),
                                                         hgtdFlexSolid, m_materialMgr->getMaterial(flexPars.material));
            GeoPhysVol* hgtdFlexPhysical = new GeoPhysVol(hgtdFlexLogical);
            flexPackagePhysical[flexVolume]->add(new GeoTransform(GeoTrf::TranslateZ3D(flexZoffset)));
            flexPackagePhysical[flexVolume]->add(hgtdFlexPhysical);
            // print out a line for each flex layer
            ATH_MSG_DEBUG( "Flex layer (" << (flexSheet ? "front" : "back") << ")" << flexSheet << ", Rmin = " << std::setw(5)
              << rInner[flexSheet] << " mm, flexZoffset = " << flexZoffset << " mm" );
            flexZoffset = flexZoffset - m_hgtdPars.flexSheetSpacing;
        }
    }

    ///////////////////
    // COOLING TUBES //
    ///////////////////

    // make list of radii of cooling tubes
    std::vector<double> coolingTubeRadii;
    double coolingTubeRadius = 130.;
    coolingTubeRadii.push_back(coolingTubeRadius);

    // two-ring layout
    if (m_geomVersion == 0) {
      ATH_MSG_INFO("Will now calculate cooling-loop positions for the two-ring layout");
      for (int i = 0; i < 18; i++) {
        coolingTubeRadius += (418-130.)/18;
        coolingTubeRadii.push_back(coolingTubeRadius);
      }
      for (int i = 0; i < 12; i++) {
        coolingTubeRadius += (658-418.)/14;
        coolingTubeRadii.push_back(coolingTubeRadius);
      }
      coolingTubeRadius = 710.;
      coolingTubeRadii.push_back(coolingTubeRadius);
      for (int i = 0; i < 7; i++) {
        coolingTubeRadius += (890-710.)/6;
        coolingTubeRadii.push_back(coolingTubeRadius);
      }
    }
    else if (m_geomVersion == 1) {
      ATH_MSG_INFO("Will now calculate cooling-loop positions for the three-ring layout");
      // inner part, even spacing from 130 mm to 674 mm, 35 rings with 16 mm spacing (first one already placed above)
      int numberOfLoops = 34;
      float loopDistance = (674.-130.)/numberOfLoops; // in mm
      for (int i = 0; i < numberOfLoops; i++) {
        coolingTubeRadius += loopDistance;
        coolingTubeRadii.push_back(coolingTubeRadius);
      }
      // outer part, even spacing from 720 mm to 900 mm, 7 rings with 30 mm spacing
      coolingTubeRadius = 720;
      coolingTubeRadii.push_back(coolingTubeRadius);
      numberOfLoops = 6;
      loopDistance = (900.-720.)/numberOfLoops;
      for (int i = 0; i < numberOfLoops; i++) {
        coolingTubeRadius += loopDistance;
        coolingTubeRadii.push_back(coolingTubeRadius);
      }
    }
    ATH_MSG_DEBUG( "Cooling tubes will be created at the following radii (" << coolingTubeRadii.size() << " in total):");
    for (size_t i = 0; i < coolingTubeRadii.size(); i++) {
        ATH_MSG_DEBUG( "   R = " << coolingTubeRadii[i] << " mm" );
    }

    ///////////////////////////////////
    // PERIPHERAL ELECTRONICS VOLUME //
    ///////////////////////////////////

    //build peripheral electronics
    GeoCylVolParams periphElPars = m_cylVolPars["HGTD::PeriphElec"];
    GeoTube*    periphElec_solid  = new GeoTube(periphElPars.rMin, periphElPars.rMax, periphElPars.zHalf);
    GeoLogVol*  periphElec_log    = new GeoLogVol(periphElPars.name, periphElec_solid, m_materialMgr->getMaterial(periphElPars.material));
    GeoPhysVol* periphElec_phys   = new GeoPhysVol(periphElec_log);

    std::array< GeoPhysVol*, 4 > moduleLayerPhysical = {}; // array of pointers to the physical volumes for the module layers which need special care

    ///////////////////////////////////////////
    // BUILD UP ALL MAIN VOLUMES IN SEQUENCE //
    ///////////////////////////////////////////

    // now build up the volumes in the order specified in the vector
    double zModuleLayerF = 0.;
    double zModuleLayerB = 0.;
    for (size_t vol = 0; vol < hgtdVolumes.size(); vol++) {

        std::string v = hgtdVolumes[vol];

        // calculate local z offsets for each main volume sequentially
        if (vol == 0) // special treatment for the first one
            m_cylVolPars[v].zOffsetLocal = motherHalfZ - m_cylVolPars[v].zHalf;

        // All but the InnerRCover, OuterRCover and peripheral cooling lines are placed relative to other components,
        // but the zOffsetLocal parameter of these volumes is left as read from the db
        else {
            if (v.substr(9,8) != "erRCover" && v != "HGTD::PeripheralCoolingLines") {
                std::string vPrev = hgtdVolumes[vol-1];
                m_cylVolPars[v].zOffsetLocal = m_cylVolPars[vPrev].zOffsetLocal - m_cylVolPars[vPrev].zHalf - m_cylVolPars[v].zHalf;
            }
        }

        // skip the tolerances - we don't actually want to create volumes for the space
        if (v.substr(0,15) == "HGTD::Tolerance") continue;

        float safety = 0.;
        if (v.substr(0,17) == "HGTD::ModuleLayer")
          safety = 10.;

        //  a disk volume to hold 4 quadrants
        GeoTube*    hgtdSubVolumeSolid    = new GeoTube(m_cylVolPars[v].rMin, m_cylVolPars[v].rMax+safety, m_cylVolPars[v].zHalf);
        GeoLogVol*  hgtdSubVolumeLogical  = new GeoLogVol(m_cylVolPars[v].name, hgtdSubVolumeSolid, m_materialMgr->getMaterial(m_cylVolPars[v].material));
        GeoPhysVol* hgtdSubVolumePhysical = new GeoPhysVol(hgtdSubVolumeLogical);

        // if building the cooling plate, also add peripheral electronics since position of those are relative to that of cooling plate
        if (v == "HGTD::CoolingPlate") {
            double zOffsetPeriphElec = m_cylVolPars[v].zHalf + periphElPars.zOffsetLocal + periphElPars.zHalf;
            // place two, one on each side of cooling plate
            for (int side = 0; side < 2; side++) {
                HGTDparent->add(new GeoTransform(GeoTrf::TranslateZ3D(m_cylVolPars[v].zOffsetLocal + pow(-1, side)*zOffsetPeriphElec)));
                HGTDparent->add(periphElec_phys);
            }

            // and the CO2 cooling tubes inside the cooling plate
            for (size_t i = 0; i < coolingTubeRadii.size(); i++) {
                // the tube itself
                GeoTorus* coolingTubeSolid = new GeoTorus(m_cylVolPars["HGTD::CoolingTubeFluid"].zHalf, m_cylVolPars["HGTD::CoolingTube"].zHalf,
                                                          coolingTubeRadii[i], 0, 2*M_PI);
                GeoLogVol* coolingTubeLogical = new GeoLogVol("HGTD::CoolingTube", coolingTubeSolid,
                                                              m_materialMgr->getMaterial(m_cylVolPars["HGTD::CoolingTube"].material));
                GeoPhysVol* coolingTubePhysical = new GeoPhysVol(coolingTubeLogical);
                hgtdSubVolumePhysical->add(coolingTubePhysical); // no transformations needed, concentric with cooling plate and centered in z
                // and the contents, i.e. the cooling fluid
                GeoTorus* coolingFluidSolid = new GeoTorus(0, m_cylVolPars["HGTD::CoolingTubeFluid"].zHalf,
                                                           coolingTubeRadii[i], 0, 2*M_PI);
                GeoLogVol* coolingFluidLogical = new GeoLogVol("HGTD::CoolingFluid", coolingFluidSolid,
                                                               m_materialMgr->getMaterial(m_cylVolPars["HGTD::CoolingTubeFluid"].material));
                GeoPhysVol* coolingFluidPhysical = new GeoPhysVol(coolingFluidLogical);
                hgtdSubVolumePhysical->add(coolingFluidPhysical); // no transformations needed, concentric with cooling plate and centered in z
            }
        }

        // module layer
        if (v.substr(0,17) == "HGTD::ModuleLayer") {

            int layer = atoi(v.substr(17,1).c_str());

            // front and back side layers are treated differently: z position of flex and module layers, and rotation
            double zFlex = 0.;
            bool Lside = layer % 2;
            if (Lside == 0) { // layers 0 and 2
                zFlex = -modulePackageHalfZtot + m_cylVolPars["HGTD::FlexPackage"].zHalf;
                zModuleLayerF = modulePackageHalfZtot - modulePackageHalfZ;
            }
            else { // layers 1 and 3
                zFlex = modulePackageHalfZtot - m_cylVolPars["HGTD::FlexPackage"].zHalf;
                zModuleLayerB = -modulePackageHalfZtot + modulePackageHalfZ;
            }

            // place flex within module packages, at different positions depending on front or back or cooling plate
            hgtdSubVolumePhysical->add(new GeoTransform(GeoTrf::TranslateZ3D(zFlex)));
            hgtdSubVolumePhysical->add(flexPackagePhysical[(Lside ? 0 : 1)]);

            float diskRotation = layer <= 1 ? m_hgtdPars.disk1Rotation : m_hgtdPars.disk2Rotation;

            HGTDparent->add(new GeoTransform( GeoTrf::TranslateZ3D(m_cylVolPars[v].zOffsetLocal) *
                                              GeoTrf::RotateZ3D(diskRotation*Gaudi::Units::deg)) );
            // one needs to check this rotation  against the "quadrot" will be used in the following

            HGTDparent->add(  hgtdSubVolumePhysical );
            moduleLayerPhysical[layer] = hgtdSubVolumePhysical;

        } // end of module package
        else {
            HGTDparent->add(new GeoTransform(GeoTrf::TranslateZ3D(m_cylVolPars[v].zOffsetLocal)));
            HGTDparent->add(hgtdSubVolumePhysical);
        }

        // print out info about each main volume
        ATH_MSG_INFO( std::setw(20) << m_cylVolPars[v].name << " ( " << std::setw(20) << m_cylVolPars[v].material
               << " ), local z = " << std::setw(6) << m_cylVolPars[v].zOffsetLocal
               << " mm, Rmin = " <<  std::setw(4) << m_cylVolPars[v].rMin
               << " mm, Rmax = " << std::setw(4) << m_cylVolPars[v].rMax
               << " mm, DZ = " << std::setw(5) << m_cylVolPars[v].zHalf << " mm" );

    } // end loop over hgtdVolumes

    ////////////////////
    // MODULE VOLUMES //
    ////////////////////

    // components for the module
    std::vector<std::string> moduleVolumes;
    moduleVolumes.push_back("HGTD::GlueAsic");
    moduleVolumes.push_back("HGTD::ASIC");
    moduleVolumes.push_back("HGTD::LGADInactive");
    moduleVolumes.push_back("SensorPlaceHolder"); // replaced below to get the numbered name right
    moduleVolumes.push_back("HGTD::GlueSensor");
    moduleVolumes.push_back("HGTD::Hybrid");
    moduleVolumes.push_back("HGTD::ModuleSpace");

    int endcap = bPos ? +2 : -2;
    double thickness = 2.*m_boxVolPars["HGTDSiSensor0"].zHalf;
    InDetDD::HGTD_ModuleDesign* moduleDesign = createHgtdDesign( thickness );

    // create the module --> each for cell and with different names
    // calculate the positions where modules should be placed in one quadrant

    int totMod = 0;
    //  this should be taken from DB or XML
    unsigned int maxRows = 21;
    if ( m_geomVersion == 0 ) maxRows = 18;

    std::array< PositionsInQuadrant, 4 > positions = prepareLayersFromQuadrants( maxRows ) ;
    // inside m_geomVersion implicitly control 3-ring layout vs 2-ring

    mirrorPositionsAroundYaxis(positions);

    for (int layer = 0; layer < 4; layer++) {
        if (m_outputIdfr) cout << "Layer #" << layer << std::endl;
        // select from front vs back side of a disk
        int Lside = layer % 2;

        std::vector<std::string> volumes = moduleVolumes;
        if ( Lside != 0 ) reverse( volumes.begin(), volumes.end() ); // reverse order of components for backside modules

        std::string sensorName = std::string("HGTDSiSensor") + std::to_string(layer);
        std::string moduleName = std::string("HGTDModule") + std::to_string(layer);

        // here we assumed all 4 layers share the same dimensions.
        // As described at HGTD_DetectorFactory::reorderRows,
        // the short edge of 4*2 cm module in the leading row is defined as local X/Width
        // in ReadoutGeometry, this short edge is also defined as Eta, since the row is roughly along radius.
        double moduleHalfWidth  = m_boxVolPars[moduleName].xHalf;    // 11m than 10 to hold wire bond
        double moduleHalfHeight = m_boxVolPars[moduleName].yHalf;

        // loop over quadrants in the current layer
        // take a prepared quadrant as protype
        PositionsInQuadrant tmpQuadrant = positions[ layer ];
        //  The relative rotation between two disks is supposed to be defined/accounted within tmpQuadrant
        for (int q = 0; q < 4; q++) {
            float quadrot = q*90.;

            for ( unsigned int row = 0; row < maxRows; row ++ ) {
                std::vector< ModulePosition > ModsPerRow = tmpQuadrant[ row ];

                // print #modules per row to fill HGTD_Identifier dictionary etc.
                if ( m_outputIdfr && q == 0 ) std::cout << " Row  #"<< row + 1 <<" :: " << ModsPerRow.size() << std::endl;

                for ( unsigned int mod = 0; mod < ModsPerRow.size(); mod ++ ) {
                    ModulePosition module = ModsPerRow[ mod ];

                    double myx = -9999999.9 , myy = -9999999.9 , myrot = -9999999.9;
                    int myphi = -1 , myeta = - 1;
                    std::string module_string = formModuleName( layer, q, maxRows, row, mod, module, myx, myy, myrot, myphi, myeta );

                    if ( module_string == "" || myrot == -9999999.9 || myeta == -1 )
                      ATH_MSG_WARNING ( " Please check the module at layer "<< layer <<" quadrant " << q <<" row "<< row <<" mod " << mod <<" not well retrieved ! " );

                    //  an hgtd module  defined in the form of  ( X, Y, Z )
                    GeoBox* moduleSolid            = new GeoBox( moduleHalfWidth, moduleHalfHeight, modulePackageHalfZ);
                    GeoLogVol* moduleLogical       = new GeoLogVol( moduleName + module_string, moduleSolid, m_materialMgr->getMaterial("std::Air"));
                    GeoFullPhysVol* modulePhysical = new GeoFullPhysVol( moduleLogical );

                    // print out one module per layer
                    if ( q == 0 && row == 0 && mod == 0 )
                      ATH_MSG_DEBUG( "Will now build up an individual HGTD module of layer " << layer << " and quadrant " << q << " (" << module_string << ")" );

                    // loop over components in module
                    for (size_t comp = 0; comp < volumes.size(); comp++) {
                        if (volumes[comp] == "SensorPlaceHolder") volumes[comp] = sensorName; // replace placeholder

                        std::string c = volumes[comp];
                        // calculate local z offsets for each sensor component sequentially
                        if (comp == 0) // special treatment for the first one
                            m_boxVolPars[c].zOffsetLocal = modulePackageHalfZ - m_boxVolPars[c].zHalf;
                        else {
                            std::string cPrev = volumes[comp-1];
                            m_boxVolPars[c].zOffsetLocal = m_boxVolPars[cPrev].zOffsetLocal - m_boxVolPars[cPrev].zHalf - m_boxVolPars[c].zHalf;
                        }

                        // skip the module space - we don't actually want to create volumes for the space
                        if (volumes[comp] == "HGTD::ModuleSpace") continue;

                        double comp_halfx = m_boxVolPars[c].xHalf;
                        double comp_halfy = m_boxVolPars[c].yHalf;

                        // and x offsets for those components that are smaller, to make room for wire bond of flex to ASIC which is larger than the sensor
                        double xOffsetLocal = moduleHalfWidth - comp_halfx;
                        // need tuning then dataBase : to make room for wire bond of flex to ASIC which is larger than the sensor

                        GeoBox* sensorCompSolidVol    = new GeoBox(comp_halfx, comp_halfy, m_boxVolPars[c].zHalf);
                        // No attachment?
                        std::string attach = (volumes[comp] == sensorName) ? "" : "_L" + std::to_string( layer ) + module_string;

                        GeoLogVol*  sensorCompLogicalVol  = new GeoLogVol( m_boxVolPars[c].name+attach, sensorCompSolidVol,
                                                                          m_materialMgr->getMaterial(m_boxVolPars[c].material));
                        GeoFullPhysVol* sensorCompPhysicalVol = new GeoFullPhysVol(sensorCompLogicalVol);

                        if (volumes[comp] == sensorName) {
                            const HGTD_ID* hgtdId = dynamic_cast<const HGTD_ID*>( m_athComps->getIdHelper() );
                            Identifier idwafer = hgtdId->wafer_id( endcap, layer, myphi, myeta );

                            // print only the first and last module of each row in the first quadrant
                            if ( q == 0 && ( mod == 0 || mod == ( ModsPerRow.size() - 1 ) ) && !m_outputIdfr ) {
                              ATH_MSG_DEBUG( "  waferHash :  " << hgtdId->wafer_hash( idwafer )
                              << " upon HGTD_ID =>  ec: " << endcap << ", layer: " << layer << ", quadrant: " << q
                              << ", row: " << myphi <<", module: "<< myeta );
                              ATH_MSG_DEBUG( " HGTD Module: " << m_boxVolPars[c].name+module_string << ", posX: " << myx << ", posY: " << myy << ", rot: " << quadrot + myrot );
                            }

                            InDetDD::HGTD_DetectorElement* detElement = new InDetDD::HGTD_DetectorElement(idwafer, moduleDesign, sensorCompPhysicalVol, m_commonItems.get());
                            m_detectorManager->addDetectorElement( detElement );

                            GeoTrf::Transform3D sensorTransform = GeoTrf::TranslateZ3D(m_boxVolPars[c].zOffsetLocal)*GeoTrf::TranslateX3D(xOffsetLocal);
                            GeoAlignableTransform* xform = new GeoAlignableTransform(sensorTransform);

                            modulePhysical->add( xform );
                            modulePhysical->add( sensorCompPhysicalVol );

                            totMod ++;
                        }
                        else {
                            modulePhysical->add(new GeoTransform(GeoTrf::TranslateZ3D(m_boxVolPars[c].zOffsetLocal)*GeoTrf::TranslateX3D(xOffsetLocal)));
                            modulePhysical->add(sensorCompPhysicalVol);
                        }

                        // print out each module component
                        if ( mod == 0 && q == 0 && volumes[comp] != sensorName )
                          ATH_MSG_DEBUG( std::setw(20) << m_boxVolPars[c].name << " ( " << std::setw(15) << m_boxVolPars[c].material
                          << " ), in-sensor-layer local z = " << std::setw(7) << m_boxVolPars[c].zOffsetLocal << " mm"
                          << ", DX = " << std::setw(5) << m_boxVolPars[c].xHalf << " mm"
                          << ", DY = " << std::setw(5) << m_boxVolPars[c].yHalf << " mm"
                          << ", DZ = " << std::setw(5) << m_boxVolPars[c].zHalf << " mm" );
                    } // end of components loop

                    double zModule = ( Lside == 0 ? zModuleLayerF : zModuleLayerB );

                    GeoTransform* moduleTransform = new GeoTransform( GeoTrf::TranslateZ3D(zModule) *
                                                                      GeoTrf::TranslateX3D(myx) *
                                                                      GeoTrf::TranslateY3D(myy) *
                                                                      GeoTrf::RotateZ3D( ( quadrot + myrot )*Gaudi::Units::deg) );
                    moduleLayerPhysical[layer]->add( moduleTransform );
                    moduleLayerPhysical[layer]->add( modulePhysical );
                } //end of modules loop
                ATH_MSG_DEBUG( "Done placing modules for row " << row );
            } // end of row loop
            ATH_MSG_DEBUG( "Done placing modules for quadrant " << q );
        } // end of quadrants loop
        ATH_MSG_DEBUG( "Done placing modules for layer " << layer );
    }

    ATH_MSG_INFO( "**************************************************" );
    ATH_MSG_INFO( "  Done building HGTD with " << totMod <<" modules " );
    ATH_MSG_INFO( "**************************************************" );

    return HGTDparent;
}

// backwards compatibility two-ring layout used for TDR studies
// careful implicit dependence on m_geomVersion to make 3-ring layout different from 2-ring
std::array< PositionsInQuadrant, 4 > HGTD_DetectorFactory::prepareLayersFromQuadrants( unsigned int nRows ) {

    std::array< PositionsInQuadrant, 4 > positions;

    // three-ring layout
    if ( m_geomVersion == 1 ) {
        // the 1st disk
        PositionsInQuadrant d0q0front = prepareQuadrantsFromRows( 0, nRows );
        reorderRows( &d0q0front );
        PositionsInQuadrant d0q0back = prepareQuadrantsFromRows( 1, nRows );
        reorderRows( &d0q0back );

        // the 2nd disk is the same as the 1st, but mirrored
        PositionsInQuadrant d1q0front = mirrorModulesInQuadrant( d0q0back );
        PositionsInQuadrant d1q0back = mirrorModulesInQuadrant( d0q0front );

        positions[ 0 ] = d0q0front; // front-side module positions
        positions[ 1 ] = d0q0back; // back-side module positions
        positions[ 2 ] = d1q0front; // front-side module positions
        positions[ 3 ] = d1q0back; // back-side module positions

    } else {
        nRows = 18; // note  21-18 = 3 elements with positions of modules in rows are left empty
        positions[0] = positions[2] = prepareQuadrantsFromRows(0, nRows ); // front-side module positions
        positions[1] = positions[3] = prepareQuadrantsFromRows(1, nRows ); // back-side module positions
    }
    return positions;
}

//  careful m_geomVersion control layout implicitly
// backward compatibility to pre-TDR two-ring layouts
// 3-ring layout differ from 2-ring here.
std::string HGTD_DetectorFactory::formModuleName( int layer, int quadrant, unsigned int maxrows, int row, int mod,
                                                  ModulePosition module,
                                                  double& myx, double& myy, double& myrot, int& phi, int& eta ) {

    std::string module_string = "";

    double x = module.x;
    double y = module.y;
    double myphi  = atan(y/x);
    double radius = std::sqrt(x*x+y*y);

    myx = radius*cos( quadrant*M_PI*0.5 + myphi );
    myy = radius*sin( quadrant*M_PI*0.5 + myphi );

    // three-ring layout
    if ( m_geomVersion == 1 ) {
        myrot = module.phiRotation;
        phi = quadrant*21 + row + 1;  // quadrant is absent ( hidden into row ) in HGTD-Identifier
        eta = mod + 1;
        //module_string = "_R" + std::to_string(phi) + "_M" + std::to_string(eta); //This was the previous string, but doesn't match expectations of HGTDSensorSD
        module_string = "_layer_" + std::to_string(layer) + "_" + std::to_string(phi) + "_" + std::to_string(eta);
    }
    // two-ring layout
    else {
        double rot = module.flipped ? 90. : 0.;
        int myrow = module.row;
        double moduleRotation = 0;
        if ( layer > 1 ) {
            myx = -myx;
            // need to rotate 180 degrees some modules in q0 and q2
            if ( quadrant%2 == 0 && myrow <= 15)  moduleRotation = 180.;
            else if (quadrant%2 == 1 && myrow > 15) moduleRotation = 180.;
        }
        eta = ( quadrant*maxrows ) + myrow;
        phi = module.el_in_row;
        myrot = moduleRotation + rot;
        module_string = "_layer_" + std::to_string(layer) + "_" + std::to_string(phi) + "_" + std::to_string(eta);
    }

    return module_string;
}

// calculate the positions of modules in a quadrant, taking care of separate schemes for two- and three-ring layouts
PositionsInQuadrant HGTD_DetectorFactory::prepareQuadrantsFromRows( int layer, unsigned int maxRow ) {

    PositionsInQuadrant rowsInQuad;
    bool isBackside = (layer % 2);
    // three-ring layout
    if ( m_geomVersion == 1 ) {
        for (size_t row = 0; row <= maxRow; row++) {
            if ( row == 13 ) continue;  // element #21 is tried since one row is skipped
            std::vector<ModulePosition> rowModulePositions = prepareModulePositionsInRowThreeRing( row, isBackside );
            rowsInQuad[ row > 13 ? row - 1 : row ] = rowModulePositions;
        }
    }
    // two-ring layout
    else {
        for (size_t row = 0; row < maxRow; row++) {
            std::vector<ModulePosition> rowModulePositions = prepareModulePositionsInRowTwoRing(row, isBackside);
            rowsInQuad[ row ] = rowModulePositions;
        }
    }
    return rowsInQuad;
}

PositionsInQuadrant HGTD_DetectorFactory::mirrorModulesInQuadrant( PositionsInQuadrant inquad ) {
    PositionsInQuadrant rowsInQuad;
    for (size_t row = 0; row < inquad.size(); row ++ ) {
        std::vector<ModulePosition> modulePositions = inquad[ row ];
        for ( size_t mod = 0; mod < modulePositions.size(); mod++ ) {
            ModulePosition old = modulePositions[mod];
            ModulePosition mirror = old;
            mirror.x = old.y;
            mirror.y = old.x;
            mirror.phiRotation = (  old.phiRotation == 0. ?  90. : 0. );
            modulePositions[mod] = mirror;
        }
        // keeping the order defined in HGTD_DetectorFactory::reorderRows
        rowsInQuad[ inquad.size() - row - 1 ] = modulePositions;
    }
    return rowsInQuad;
}

// calculate module positions for the three-ring layout, updated to agree with post-TDR developments
std::vector< ModulePosition > HGTD_DetectorFactory::prepareModulePositionsInRowThreeRing( int row, int back ) {

  // below parameters should be collected into xml or RDB
  int index_XYcoord_change = 14;

  // height is the short edge of module, width is the longer edge
  float halfWidth = .5*40., halfHeight = .5*21.8; // bare module
  float midR = 230., midR2 = 470.5, maxRcut = 660., maxOuterR = 670.;
  // special tweak needed to avoid overlapping holes on front and back
  if (row == 21 and back==1) {
    midR2 = 510.;
  }
  float readoutRowSpace = 1.0;
  bool extrude = ( ( row == 6 || row == 18 ) && !back  ) || // front side
                 ( ( row == 2 || row == 11 || row == 12 || row == 17 ) && back ); // back side

  // instead of attempting to re-calculate the leading module per row, just pick up from dataBase,
  // numbers here taken from spreadsheet at https://cernbox.cern.ch/index.php/s/PPXEWSBnBjwI7UU
  std::array< float, 22 > ModStarting = { 122., 122.7, 89.85, 123.5, 175.4, 257.4, 287.5, 298.4, 287.5, 304.5, 287.5, 304.5, 287.5, 0.0, 299.7,
                                          130., 114.7, 131.45, 164.45, 216.35, 205.45, 257.35 };

  std::array< float, 22 > ModStartBack = { 130., 114.7, 97.85, 131.5, 164.5, 246.5, 298.4, 287.5, 298.4, 287.5, 304.5, 287.5, 304.5, 0.0, 287.5,
                                           122., 122.7, 123.45, 172.45, 205.45, 216.35, 246.45 };
  /*  row == 13 will be skipped from outside, and, since row == 15 XY flip take place. */

  // this is a bit hacky for this layout, a db-based solution will help in the future
  int useCorner = 0;
  if ( (  ( row == 1 || row == 5 || row == 15 || row == 19 ) && ! back ) || // front side
       (  ( row == 0 || row == 8 || row == 16 || row == 18 || row == 21 ) && back ) // back side
     )  useCorner = 1;
  if ( row == 17 ) useCorner = 2;
  // in some exceptional cases the spacing will be smaller even though the module crossed the ring boundary
  float backshift = 6.;

  // the new layout tune makes small adjustments (usually 2~3 mm) for the last modules of some rows.
  // even though most of element is zero for most of modules, we store these adjustments in a 2D array for now
  float tailModCorrection[ 22 ][ 19 ];
  for ( int r = 0; r < 22; r ++ )
    for ( int m = 0; m < 19; m ++ ) tailModCorrection[r][m] = 0.;
  tailModCorrection[11][4] = tailModCorrection[12][2] = 10.;

  // TDR layout: spaceSmallR = 5.5 , spaceMediumR = 8.4 , spaceLargeR = 14.5
  float spaceSmallR = 3.7 , spaceMediumR = 6.6, spaceLargeR = 12.7; // updated spacings from post-TDR developments

  float backsideSmallR = spaceSmallR;
  float backsideMediumR = spaceMediumR;
  float backsideLargeR = spaceLargeR;

  float extendedWidth = readoutRowSpace + 2.*halfWidth;

  float posRadius = 0.;
  float posOfLastPlacedModule = 0.; // start with dummy position
  int moduleCounter = 0;
  std::vector< ModulePosition > rowModulePositions;

  float effectiveRow = row;
  // note the flipping of effectiveRow, it is the cause of XY flipping in later occurrance
  if ( row == index_XYcoord_change ) effectiveRow = 13;
  if ( row > index_XYcoord_change ) effectiveRow -= ( index_XYcoord_change + 1 );

  // x coordinate for vertical rows
  float rowCentPos = 0.5*extendedWidth*( 2*effectiveRow + 1 );

  if ( extrude )  maxRcut = maxOuterR;
  while ( posRadius < maxRcut ) {
    // horizontal rows need care (restart from other edge of quadrant), this variable helps get the inner radius right
    // y coordinate for vertical rows, then x coordinate for modules
    float modPos_row = -999.;

    // for the first module in each row, take the starting position from the arrays created earlier from the spreadsheet
    if ( moduleCounter == 0 )  { // leading module per row
      modPos_row = ( back ?  ModStartBack[row] : ModStarting[row] );
      modPos_row += halfHeight;
    }
    // the rest of the modules follow sequential, radius-dependent placement rules
    else {
      float prevX = rowModulePositions[ moduleCounter - 1 ].x;
      float prevY = rowModulePositions[ moduleCounter - 1 ].y;
      float spacing = back ? backsideSmallR : spaceSmallR;

      // increase the spacing by the ring it will fallin
      float ringCrossRcorner =  std::sqrt( ( prevY + halfHeight)*( prevY + halfHeight) +
                       ( prevX + halfWidth )*( prevX + halfWidth ) );
      float ringCrossRcenter =  std::sqrt( prevY*prevY + prevX*prevX );

      bool tuned_center = ( row == 3 && (  moduleCounter == 3 && !back ) ) || // front, row 3
                      ( row == 20 && moduleCounter == 8 && !back ) || // front, row 20
                      ( row == 21 && moduleCounter == 6 && back ); // back, row 21
      if ( useCorner == 2 ) {
        if ( ( moduleCounter == 3 && ! back ) || ( ( moduleCounter ==  3 || moduleCounter == 4 ) && back ) ) {
          ringCrossRcenter -= backshift;
          if ( ringCrossRcenter > midR && ringCrossRcenter <= midR2 ) spacing = back ? backsideMediumR : spaceMediumR;
          if ( ringCrossRcenter > midR2 ) spacing = back ? backsideLargeR : spaceLargeR;
        }
    else {
          if ( ringCrossRcorner > midR && ringCrossRcorner <= midR2 ) spacing = back ? backsideMediumR : spaceMediumR;
          if ( ringCrossRcorner > midR2 ) spacing = back ? backsideLargeR : spaceLargeR;
        }
      }
      else if ( useCorner == 1 ) {
        if ( ringCrossRcorner > midR && ringCrossRcorner <= midR2 ) spacing = back ? backsideMediumR : spaceMediumR;
        if ( ringCrossRcorner > midR2 ) spacing = back ? backsideLargeR : spaceLargeR;
      }
      else {
        if ( tuned_center ) ringCrossRcenter -= backshift;
        if ( ringCrossRcenter > midR && ringCrossRcenter <= midR2 ) spacing = back ? backsideMediumR : spaceMediumR;
        if ( ringCrossRcenter > midR2 ) spacing = back ? backsideLargeR : spaceLargeR;
      }

      modPos_row = posOfLastPlacedModule + 2.*halfHeight + spacing;

      if ( back ) modPos_row -= tailModCorrection[ row ][ moduleCounter ];

    } // endif  non-leading module

    // check and limit the length of the row
    posRadius = std::sqrt( ( rowCentPos + halfWidth )*( rowCentPos + halfWidth ) +
                           ( modPos_row + halfHeight)*( modPos_row + halfHeight ) );
    if ( posRadius  > maxRcut ) {
      ATH_MSG_DEBUG(" row " << row <<" finished with " << moduleCounter <<"  modules ");
      break;
    }

    // the X and Y coordinates need to be flipped if this row is horizontal,
    // needed only for backwards compatibility for two-ring layout
    ModulePosition modu =  { modPos_row, rowCentPos, 0., false, row, moduleCounter };
    ModulePosition moduFlipped =  { rowCentPos, modPos_row, 90., true, row, moduleCounter };

    // eventually arrived at a module to append into this row
    if ( row > index_XYcoord_change ) rowModulePositions.push_back( modu );
    else rowModulePositions.push_back( moduFlipped );

    // the spreadsheet gave the center of bottom edge of a module, so an adjustment by halfHeight is needed
    if ( m_outputIdfr ) ATH_MSG_DEBUG( " Row " << ( row <= index_XYcoord_change ? effectiveRow + 1 : 36 - row )
      << " Module " << moduleCounter + 1 <<" at (x,y) : "
      << ( row > index_XYcoord_change ? rowModulePositions.back().x - halfHeight : rowModulePositions.back().x ) << ", "
      << ( row > index_XYcoord_change ? rowModulePositions.back().y : rowModulePositions.back().y - halfHeight ) );

    posOfLastPlacedModule = modPos_row;
    moduleCounter ++;
  }  // endof while loop

  if ( m_outputIdfr ) std::cout << "Total #Module " << rowModulePositions.size() <<" at original row " << row << std::endl;

  return rowModulePositions;
}

// adjust row ordering to adhere to a convention in upper right quadrant of front :
// 0 is first/bottom of many horizontal laid module along X-coord ( Y is almost zero )
// increasing as they're stacked upwards ( anticlockwise ), the upmost horizontal is 6
// then the outermost ( right-most ) vertical is 7, increase as they are laid to left
// 20 is the last ( left ) adhereing to Y-coord ( X is almost zero )
//
int HGTD_DetectorFactory::reorderRows( PositionsInQuadrant* quadrant ) {

    PositionsInQuadrant tmpQuadrant;
    int xchng = 0;
    unsigned int numrow =  quadrant->size();
    for ( unsigned int r = 0; r < numrow; r ++ ) {
        unsigned int idx = r > 13 ?  13 + numrow - r : r;
        ATH_MSG_DEBUG( " original row " << ( r <= 12 ? r : r + 1 ) <<" new row " << idx + 1
                        << " : "<< numrow  );
        tmpQuadrant[ idx ] = quadrant->at( r );
        if ( idx != r ) xchng++;
    }
    // refill the quadrant
    for ( unsigned int r = 0; r < numrow; r++ )
    quadrant->at( r ) = tmpQuadrant[ numrow - 1 - r ];

    return xchng;
}

void HGTD_DetectorFactory::mirrorPositionsAroundYaxis(std::array< PositionsInQuadrant, 4 >& arr) {
  for (auto& layer : arr) {
    for (auto& row : layer) {
      for (auto& module : row) {
        module.x = -module.x;
      }
    }
  }
}


InDetDD::HGTD_ModuleDesign* HGTD_DetectorFactory::createHgtdDesign( double thickness ) {

    double phiPitch = 1.3; // mm
    double etaPitch = 1.3; // mm

    // for each module there are two ASICs connected, need further confirmation for below specifications
    // the "collumns" along module short edge align to a direction just following the EndOfCollumn(EoC) of each ASIC.
    int circuitsPerColumn = 1;
    int circuitsPerRow = 2;  // along the long edge of module, tentatively following IBL descriptin

    // each ASIC just cover half ( 15*15 pixels ) of a module
    int cellColumnsPerCircuit = 15;
    int cellRowsPerCircuit = 15;   //  along module FULL length there are 2*15  rows

    int diodeColumnsPerCircuit = cellColumnsPerCircuit;
    int diodeRowsPerCircuit = cellRowsPerCircuit;

    std::shared_ptr<const PixelDiodeMatrix> normalCell = InDetDD::PixelDiodeMatrix::construct(phiPitch, etaPitch);
    std::shared_ptr<const PixelDiodeMatrix> singleRow  = InDetDD::PixelDiodeMatrix::construct(InDetDD::PixelDiodeMatrix::phiDir, 0,
                                                                                              normalCell, diodeColumnsPerCircuit, 0);
    std::shared_ptr<const PixelDiodeMatrix> fullMatrix = InDetDD::PixelDiodeMatrix::construct(InDetDD::PixelDiodeMatrix::etaDir, 0,
                                                                                              singleRow, 2*diodeRowsPerCircuit, 0); // note 30 = 2*15 rows adopted

    DetectorDesign::Axis yDirection = InDetDD::DetectorDesign::yAxis;

    InDetDD::HGTD_ModuleDesign* design = new InDetDD::HGTD_ModuleDesign(thickness,
                                                                        circuitsPerColumn, circuitsPerRow,
                                                                        cellColumnsPerCircuit, cellRowsPerCircuit,
                                                                        diodeColumnsPerCircuit, diodeRowsPerCircuit,
                                                                        fullMatrix,
                                                                        InDetDD::CarrierType::electrons, 1, yDirection );

    return design;
}

const HGTD_DetectorManager* HGTD_DetectorFactory::getDetectorManager() const {
    return m_detectorManager;
}

// backwards compatibility two-ring layout used for TDR studies
std::vector<ModulePosition> HGTD_DetectorFactory::prepareModulePositionsInRowTwoRing(int row, bool back) {

    unsigned int module = 0;
    std::vector<ModulePosition> modulePositions;
    double posOfLastPlacedModule = 0.;

    while (true) {
        //ATH_MSG_INFO) << "Will now place module " << module );
        // horizontal rows need care (restart from other edge of quadrant), this variable helps get the inner radius right
        // in quadrant 0 - ie top right quadrant
        // row 0 = bottom horizontal row. numbering grows upwards and counterclockwise; row 17=leftmost vertical row
        // rowForInnerRadius = 0-1 for vertical rows too
        int rowForInnerRadius = row; // because row 16-17 are equivalent to 0-1 regarding module placement
        if (row == 17) rowForInnerRadius = 0;
        if (row == 16) rowForInnerRadius = 1;

        // params needed frequently below
        double moduleWidth = m_boxVolPars["HGTDModule0"].yHalf*2;
        double moduleHeight = m_boxVolPars["HGTDModule0"].xHalf*2;
        double rInner = m_cylVolPars["HGTD::ModuleLayer0"].rMin;
        double rMid = m_hgtdPars.rMid;
        double rOuter = m_hgtdPars.rOuter;
        double rowSpaceSide = m_hgtdPars.rowSpaceSide;

        // x coordinate for vertical rows
        double rowCenterPos = (moduleWidth + rowSpaceSide)*(rowForInnerRadius + 0.5);
        // y coordinate for vertical rows
        double modulePosAlongRow = -99.; // mock value for now
        // for the first module, pick the right starting point
        if (module == 0) {
            // start at inner radius and include any offset for backside.
            if (rowForInnerRadius < 3) {
            modulePosAlongRow = std::sqrt( pow(rInner, 2) - pow((moduleWidth + rowSpaceSide)*rowForInnerRadius, 2) )
                                  + back*m_hgtdPars.rowBacksideInnerShift + moduleHeight/2;
            }
            else { // later modules start at 2*moduleWidth, with offset for the backside
                double backSpacing = m_hgtdPars.rowBacksideInnerShift;
                if (back && (rowCenterPos - moduleWidth/2 > rMid)) {
                    backSpacing = m_hgtdPars.rowBacksideOuterShift;
                }
                modulePosAlongRow = 2*(moduleWidth + rowSpaceSide) + moduleHeight/2 + back*backSpacing;
            }
        }
        // the rest of the modules follow sequential, radius-dependent placement rules
        else {
            // cppcheck-suppress containerOutOfBounds; false positive
            ModulePosition prev = modulePositions.back();
            double spacing = m_hgtdPars.moduleSpaceInner;
            // if the previous module was completely outside rMid, increase the spacing
            // (+1 mm is a needed shift for full coverage - could need tweaking if layout parameters change!)
            float innermostCornerR = std::sqrt( pow(prev.y - moduleHeight/2, 2) + pow(prev.x - moduleWidth/2, 2) ) + 1.;
            if (innermostCornerR > rMid) {
                spacing = m_hgtdPars.moduleSpaceOuter;
            }
            // for the back the large spacing starts as soon as the space would entirely be outside R = 320 mm
            if (back) {
            double startOfSpaceAlongRow = std::sqrt( pow(prev.y + moduleHeight/2, 2) + pow(prev.x - moduleWidth/2, 2) ) - 2;
                if (startOfSpaceAlongRow > rMid) {
                    spacing = m_hgtdPars.moduleSpaceOuter;
                }
            }
            // correction to the first two spaces on the front side, compensating for the 2 mm introduced at beginning of a row
            if (!back && rowForInnerRadius < 8 && module < 3) {
                spacing -= 1.; // should result in 3 mm instead of 4 mm for innermost spaces on rows starting at R = 120
            }
            // squeeze in a few more modules at the end of some rows
            double maxRcut = rOuter+20;
            if(row == 8 || row == 9 || row == 10) {
                maxRcut = 661;
                if(row == 8 && module > 12) spacing -= 4;
            }
            else if ( row == 11 && module > 9) {
                maxRcut = 662;
                spacing -= 6;
            }
            else if (row == 12 && back) {
                maxRcut = 665;
            }
            else if(row == 13 && module > 5) {
                maxRcut = 666;
                if (!back && module > 6 ) spacing -= 8.5;
                else if (back && module > 5) spacing -= 2;
            }
            else if (row == 14 && module > 3) {
                maxRcut = 665;
                spacing -= 5;
            }
            else if (row == 15) {
                maxRcut = 669;
                spacing -= 5.5;
            }
            modulePosAlongRow = posOfLastPlacedModule + moduleHeight + spacing;
            // stop if the next module will extend outside the max allowed radius
            // HC max radius is 665
            if ( std::sqrt( pow(rowCenterPos + moduleWidth/2, 2) + pow(modulePosAlongRow + moduleHeight/2, 2) ) > maxRcut) {
                break;
            }
        }
        ModulePosition m = {rowCenterPos, modulePosAlongRow, 0, true, row, (int)module};
        modulePositions.push_back(m);
        posOfLastPlacedModule = modulePosAlongRow;
        module += 1;
    } // end of loop over modules

    // finally, flip x and y for all modules if this row is horizontal
    if (row < 16) {
        // ATH_MSG_INFO) << "Flipping x and y for modules in row " << row );
        for (size_t i=0; i < modulePositions.size(); i++) {
            ModulePosition old = modulePositions[i];
            ModulePosition rotated = old;
            rotated.x = old.y;
            rotated.y = old.x;
            rotated.flipped = !old.flipped;
            rotated.row = old.row;
            rotated.el_in_row = old.el_in_row;
            modulePositions[i] = rotated;
        }
    }

    ATH_MSG_DEBUG( "row = " << row );
    for(size_t i=0; i < modulePositions.size(); i++) {
        ATH_MSG_DEBUG( "Module " << i << " at (x,y) = (" << modulePositions[i].x << "," << modulePositions[i].y << ")" );
    }

    return modulePositions;
}
