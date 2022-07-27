/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "OraclePixGeoManager.h"
#include "PixelLegacyManager.h"
#include "PixelMaterialMap.h"
#include "PixelStaveTypes.h"
#include "InDetGeoModelUtils/PairIndexMap.h"
#include "InDetGeoModelUtils/TopLevelPlacements.h"
#include "InDetGeoModelUtils/InDetMaterialManager.h"

// to permit access to StoreGate
#include "StoreGate/StoreGateSvc.h"

#include "GeometryDBSvc/IGeometryDBSvc.h"
#include "GeoModelInterfaces/IGeoDbTagSvc.h"
#include "GeoModelUtilities/DecodeVersionKey.h"
#include "GeoModelKernel/GeoMaterial.h"
#include "GeoModelKernel/Units.h"

//
// Get the pixelDD Manager from SG.
//
#include "PixelReadoutGeometry/PixelDetectorManager.h"
#include "RDBAccessSvc/IRDBAccessSvc.h"
#include "RDBAccessSvc/IRDBRecord.h"
#include "RDBAccessSvc/IRDBRecordset.h"

//
// Distorted material manager
//
#include "InDetGeoModelUtils/DistortedMaterialManager.h"
#include "GaudiKernel/SystemOfUnits.h"

using InDetDD::PixelDetectorManager; 

#include <string>
#include <iostream>
#include <sstream>

using namespace std;

OraclePixGeoManager::OraclePixGeoManager(PixelGeoModelAthenaComps * athenaComps)
  : PixelGeometryManager(athenaComps), 
    m_eta(0),
    m_phi(0),
    m_currentLD(0),
    m_BarrelEndcap(0),
    m_side(0),
    m_diskFrontBack(0),
    m_servicesOnLadder(true),  
    m_services(true),  
    m_initialLayout(false), 
    m_dc1Geometry(false),
    m_alignable(true),
    m_ibl(false),
    m_PlanarModuleNumber(0),
    m_3DModuleNumber(0),
    m_dbm(false),
    m_legacyManager(nullptr),
    m_gangedIndexMap(nullptr),
    m_frameElementMap(nullptr),
    m_diskRingIndexMap(nullptr),
    m_zPositionMap(nullptr),
    m_dbVersion(0),
    m_defaultLengthUnit(Gaudi::Units::mm)
{
  m_commonItems = nullptr;
  m_pDDmgr = nullptr;

  init();
}

void
OraclePixGeoManager::init()
{
  if (msgLvl(MSG::DEBUG)) msg(MSG::DEBUG) << "Using ORACLE PIXEL GEOMETRY MANAGER" << endmsg;

  IRDBAccessSvc *rdbSvc = athenaComps()->rdbAccessSvc();
  const IGeoDbTagSvc *geoDbTag = athenaComps()->geoDbTagSvc();

  // Get version tag and node for Pixel.
  DecodeVersionKey versionKey(geoDbTag,"Pixel");
  const std::string& detectorKey  = versionKey.tag();
  const std::string& detectorNode = versionKey.node();

  // Get version tag and node for InnerDetector.
  DecodeVersionKey indetVersionKey(geoDbTag,"InnerDetector");

  m_versionTag = rdbSvc->getChildTag("Pixel", versionKey.tag(), versionKey.node());

/////////////////////////////////////////////////////////
//
// Gets the structures from the det store
//
/////////////////////////////////////////////////////////

  if(msgLvl(MSG::INFO)) {
    msg(MSG::INFO) << "Retrieving Record Sets from database ..." << endmsg;
    msg(MSG::INFO) << "Key = " << detectorKey << " Node = " << detectorNode << endmsg;
  }

  //atls = rdbSvc->getRecordset("AtlasMother",geoModel->atlasVersion(), "ATLAS");
  m_PixelSwitches      = rdbSvc->getRecordsetPtr("PixelSwitches",          detectorKey, detectorNode);
  m_PixelBarrelGeneral = rdbSvc->getRecordsetPtr("PixelBarrelGeneral",     detectorKey, detectorNode);
  m_PixelBarrelService = rdbSvc->getRecordsetPtr("PixelBarrelService",     detectorKey, detectorNode);
  m_PixelCommon        = rdbSvc->getRecordsetPtr("PixelCommon",            detectorKey, detectorNode);
  m_PixelEnvelope      = rdbSvc->getRecordsetPtr("PixelEnvelope",          detectorKey, detectorNode);
  m_PixelDisk          = rdbSvc->getRecordsetPtr("PixelDisk",              detectorKey, detectorNode);
  m_PixelDiskRing      = rdbSvc->getRecordsetPtr("PixelDiskRing",          detectorKey, detectorNode);
  m_PixelRing          = rdbSvc->getRecordsetPtr("PixelRing",              detectorKey, detectorNode);
  m_PixelEndcapGeneral = rdbSvc->getRecordsetPtr("PixelEndcapGeneral",     detectorKey, detectorNode);
  m_PixelEndcapService = rdbSvc->getRecordsetPtr("PixelEndcapService",     detectorKey, detectorNode);
  m_PixelEnvelopeService = rdbSvc->getRecordsetPtr("PixelEnvelopeService",     detectorKey, detectorNode);
  m_PixelLayer         = rdbSvc->getRecordsetPtr("PixelLayer",             detectorKey, detectorNode);
  m_PixelModule        = rdbSvc->getRecordsetPtr("PixelModule",            detectorKey, detectorNode);
  m_PixelModuleSvc     = rdbSvc->getRecordsetPtr("PixelModuleSvc",            detectorKey, detectorNode);
  m_PixelStave         = rdbSvc->getRecordsetPtr("PixelStave",             detectorKey, detectorNode);
  m_PixelStaveZ        = rdbSvc->getRecordsetPtr("PixelStaveZ",            detectorKey, detectorNode);
  m_PixelTopLevel      = rdbSvc->getRecordsetPtr("PixelTopLevel",          detectorKey, detectorNode);
  m_PixelReadout       = rdbSvc->getRecordsetPtr("PixelReadout",           detectorKey, detectorNode);
  m_PixelGangedPixels  = rdbSvc->getRecordsetPtr("GangedPixels",           detectorKey, detectorNode);
  m_PixelBarrelCable   = rdbSvc->getRecordsetPtr("PixelBarrelCable",       detectorKey, detectorNode);
  m_PixelTMT           = rdbSvc->getRecordsetPtr("PixelTMT",               detectorKey, detectorNode);
  m_PixelOmega         = rdbSvc->getRecordsetPtr("PixelOmega",             detectorKey, detectorNode);
  m_PixelOmegaGlue     = rdbSvc->getRecordsetPtr("PixelOmegaGlue",         detectorKey, detectorNode);
  m_PixelAlTube        = rdbSvc->getRecordsetPtr("PixelAlTube",            detectorKey, detectorNode);
  m_PixelFluid         = rdbSvc->getRecordsetPtr("PixelFluid",             detectorKey, detectorNode);
  m_PixelConnector     = rdbSvc->getRecordsetPtr("PixelConnector",         detectorKey, detectorNode);
  m_PixelPigtail       = rdbSvc->getRecordsetPtr("PixelPigtail",           detectorKey, detectorNode);
  m_PixelSimpleService = rdbSvc->getRecordsetPtr("PixelSimpleService",     detectorKey, detectorNode);
  m_PixelFrame         = rdbSvc->getRecordsetPtr("PixelFrame",             detectorKey, detectorNode);
  m_PixelFrameSect     = rdbSvc->getRecordsetPtr("PixelFrameSect",         detectorKey, detectorNode);
  m_PixelIBLStave      = rdbSvc->getRecordsetPtr("PixelIBLStave"  ,        detectorKey, detectorNode);
  m_PixelIBLSupport    = rdbSvc->getRecordsetPtr("PixelIBLSupport",        detectorKey, detectorNode);
  m_PixelIBLFlex    = rdbSvc->getRecordsetPtr("PixelIBLFlex",        detectorKey, detectorNode);
  m_PixelIBLFlexMaterial    = rdbSvc->getRecordsetPtr("PixelIBLFlexMaterial",        detectorKey, detectorNode);
  m_PixelIBLGlueGrease    = rdbSvc->getRecordsetPtr("PixelIBLGlueGrease",        detectorKey, detectorNode);
  m_PixelConicalStave     = rdbSvc->getRecordsetPtr("PixelConicalStave",             detectorKey, detectorNode);

  // Weights table 
  m_weightTable = rdbSvc->getRecordsetPtr("PixelWeights", detectorKey, detectorNode);

  // Extra Scaling table. This is used for extra material studies. For nominal material the table should be empty.
  // NB this is at InnerDetector level node.
  m_scalingTable = rdbSvc->getRecordsetPtr("PixelMatScaling", indetVersionKey.tag(), indetVersionKey.node());

  // MaterialMap
  m_materialTable = rdbSvc->getRecordsetPtr("PixelMaterialMap", detectorKey, detectorNode);

  // Pixel stave types
  m_staveTypeTable = rdbSvc->getRecordsetPtr("PixelStaveType", detectorKey, detectorNode);

  // DBM
  m_DBMTelescope  = rdbSvc->getRecordsetPtr("DBMTelescope",      detectorKey, detectorNode);
  m_DBMBracket    = rdbSvc->getRecordsetPtr("DBMBracket",        detectorKey, detectorNode);
  m_DBMCage    = rdbSvc->getRecordsetPtr("DBMCage",        detectorKey, detectorNode);
  m_DBMModule     = rdbSvc->getRecordsetPtr("DBMModule",         detectorKey, detectorNode);
  m_dbmWeightTable  = rdbSvc->getRecordsetPtr("DBMMaterials",      detectorKey, detectorNode);

  m_dbVersion = determineDbVersion();

  if(msgLvl(MSG::DEBUG)) msg(MSG::DEBUG) << "Database version number: " << m_dbVersion << endmsg;

  if (m_dbVersion < 4) {
    m_legacyManager = new PixelLegacyManager(rdbSvc,  detectorKey, detectorNode);
  }

  if(msgLvl(MSG::INFO)) msg(MSG::INFO) << "... Record Sets retrieved." << endmsg;

  m_distortedMatManager = new InDetDD::DistortedMaterialManager;
 
  // Set default lenth unit to Gaudi::Units::mm for newer version and Gaudi::Units::cm for older versions
  m_defaultLengthUnit =  (m_dbVersion < 3) ? Gaudi::Units::cm : Gaudi::Units::mm;

  // Get the top level placements
  m_placements = new TopLevelPlacements(m_PixelTopLevel);

  // If all individual pieces are not present, then actually all are present.
  m_allPartsPresent = (!m_placements->present("Barrel") && !m_placements->present("EndcapA") &&  !m_placements->present("EndcapC"));

  // cache the number of inner frames
  if (m_dbVersion < 3) {
    m_barrelInFrames =  (*m_PixelBarrelGeneral)[0]->getInt("NFRAMEIN");
    m_endcapInFrames =  (*m_PixelEndcapGeneral)[0]->getInt("NFRAMEIN");
  } else {
    m_barrelInFrames = 0;
    m_endcapInFrames = 0;
  }

  //
  // Get the InDet material manager. This is a wrapper around the geomodel one with some extra functionality to deal
  // with weights table if it exists
 
  m_pMatMgr = new InDetMaterialManager("PixelMaterialManager", athenaComps());
  m_pMatMgr->addWeightTable(m_weightTable, "pix");
  m_pMatMgr->addScalingTable(m_scalingTable);

  // add the DBM weight table
  m_pMatMgr->addWeightTable(m_dbmWeightTable, "pix");

  // Create material map
  m_materialMap = new PixelMaterialMap(db(), m_materialTable);
  if (m_materialTable->size() == 0) addDefaultMaterials();

  // Create stave type map
  m_pixelStaveTypes = new PixelStaveTypes(db(), m_staveTypeTable);
 
  
  //
  // Print the version number for the barrel and endcap geometry
  //  
  //cout << "Instantiating Pixel Detector" << endl;
  //cout << "Barrel Version " << this->PixelBarrelMajorVersion() << "." << this->PixelBarrelMinorVersion() << endl;
  //cout << "Endcap Version " << this->PixelEndcapMajorVersion() << "." << this->PixelEndcapMinorVersion() << endl;
}

InDetMaterialManager* OraclePixGeoManager::getMaterialManager()
{
  return m_pMatMgr;
}

PixelLegacyManager * OraclePixGeoManager::legacyManager()
{
  return m_legacyManager;
}


OraclePixGeoManager::~OraclePixGeoManager()
{
  delete m_placements;
  delete m_distortedMatManager;
  delete m_materialMap;
  delete m_pixelStaveTypes;
  delete m_legacyManager;
  delete m_pMatMgr;
  delete m_gangedIndexMap;
  delete m_diskRingIndexMap;
  delete m_zPositionMap;
  delete m_frameElementMap;
}


InDetDD::SiCommonItems * 
OraclePixGeoManager::commonItems()
{
  return m_commonItems;
}


const InDetDD::SiCommonItems * 
OraclePixGeoManager::commonItems() const
{
  return m_commonItems;
}


void 
OraclePixGeoManager::setCommonItems(InDetDD::SiCommonItems * commonItems)
{
  m_commonItems = commonItems;
}


const PixelID * 
OraclePixGeoManager::getIdHelper() 
{
  return athenaComps()->getIdHelper();
}


const GeoTrf::Transform3D & 
OraclePixGeoManager::partTransform(const std::string & partName) const 
{
  return m_placements->transform(partName);
}


bool 
OraclePixGeoManager::partPresent(const std::string & partName) const
{
  // First check if overridden from text file.
  if (partName == "Barrel") {
    if (db()->testField("PixelCommon","DOBARREL")) {
      return db()->getInt("PixelCommon","DOBARREL");
    }
  } else if (partName == "EndcapA" || partName == "EndcapC") {
    if (db()->testField("PixelCommon","DOENDCAPS")) {
      return db()->getInt("PixelCommon","DOENDCAPS");
    }
  }
  // otherwise check database.
  return (m_allPartsPresent || m_placements->present(partName));
}

/////////////////////////////////////////////////////////
//
// Setting of Layer/Disk and Barrel/EndCap
//
/////////////////////////////////////////////////////////
void OraclePixGeoManager::SetCurrentLD(int i)
{
  if(isBarrel() ) {
    if(i <= PixelBarrelNLayer()) {
      m_currentLD=i;
    } else {
      msg(MSG::ERROR) << "Layer set out of bounds: " << i << ", Setting it to 0" << endmsg;
      m_currentLD = 0;
    }
    if (msgLvl(MSG::DEBUG)) msg(MSG::DEBUG) <<" Current layer set to " << m_currentLD << endmsg;
  }
  if(isEndcap() ) {
     if (i<= PixelEndcapNDisk() ) {
       m_currentLD=i;
     } else {
       msg(MSG::ERROR) << "Disk set out of bounds: " << i << ", Setting it to 0" << endmsg;
       m_currentLD = 0;
     }
     if (msgLvl(MSG::DEBUG)) msg(MSG::DEBUG) <<" Current disk set to " << m_currentLD << endmsg;
  } 
 if(isDBM() ) {
     if (i<= 2 ) {
       m_currentLD=i;
     } else {
       msg(MSG::ERROR) << "DBM: Disk set out of bounds: " << i << ", Setting it to 0" << endmsg;
       m_currentLD = 0;
     }
     if (msgLvl(MSG::DEBUG)) msg(MSG::DEBUG) <<" Current disk set to " << m_currentLD << endmsg;
  }
}

void OraclePixGeoManager::SetBarrel() {
  //msg(MSG::DEBUG) << "Setting Barrel" << endmsg;
  m_BarrelEndcap = 0;
}
void OraclePixGeoManager::SetEndcap() {
  m_BarrelEndcap = 1;
  //msg(MSG::DEBUG) << "Setting Endcap" << endmsg;
}
void OraclePixGeoManager::SetPartsDBM() {
  m_BarrelEndcap = 2;
  //msg(MSG::DEBUG) << "Setting DBM" << endmsg;
}
/////////////////////////////////////////////////////////
//
// Check if the current layer/disk has to be retrieved
// By default in MySQL all the three layers have to be retrieved
// If m_initialLayout is true we ignore NOVA for layer/disk 1.
// For initial layout layer=1 and disk=1 (2nd layer and disk) is missing.
//
/////////////////////////////////////////////////////////
bool OraclePixGeoManager::isLDPresent() {
  if(isBarrel()) {
    if (m_initialLayout && m_currentLD == 1) return false;
    std::ostringstream A;
    A << "_" << m_currentLD;
    // More than 3 layers not yet supported in database so
    // if not present in text file assume using this layer
    return db()->getInt(m_PixelBarrelGeneral,"USELAYER"+A.str());
  }
  if(isEndcap() ) {
    if (m_initialLayout && m_currentLD == 1) return false;
    std::ostringstream A;
    A << "_" << m_currentLD;
    // More than 3 disks not yet supported in database so
    // if not present in text file assume using this disks
    return db()->getInt(m_PixelEndcapGeneral,"USEDISK"+A.str());
  }
  return false;
}


bool OraclePixGeoManager::isBarrel() {
  return m_BarrelEndcap == 0;
}
bool OraclePixGeoManager::isEndcap() {
  return m_BarrelEndcap == 1;
  return false;
}
bool OraclePixGeoManager::isDBM() {
  return m_BarrelEndcap == 2;
}

bool OraclePixGeoManager::DoServices() {
  return m_services;
}
bool OraclePixGeoManager::DoServicesOnLadder() {
  return m_servicesOnLadder;
}

bool OraclePixGeoManager::InitialLayout() const {
  return m_initialLayout;
}

void OraclePixGeoManager::SetDC1Geometry(bool flag) {
  m_dc1Geometry = flag; 
  if (m_legacyManager) m_legacyManager->SetDC1Geometry(flag);
}

bool OraclePixGeoManager::DC1Geometry() const {
  return m_dc1Geometry;
}

bool OraclePixGeoManager::Alignable() const {
  return m_alignable;
}


PixelDetectorManager* OraclePixGeoManager::GetPixelDDManager() {
  if(m_pDDmgr == nullptr) {
    //
    // retrieve the pointer to the DD manager
    //
    StatusCode sc = athenaComps()->detStore()->retrieve(m_pDDmgr);  
    if (sc.isFailure()) {
      msg(MSG::ERROR) << "Cannot retrieve PixelDetectorManager" << endmsg;
    } 
  }
  return m_pDDmgr;
}  


InDetDD::DistortedMaterialManager *
OraclePixGeoManager::distortedMatManager() {
  return m_distortedMatManager;
}  


/////////////////////////////////////////////////////////
//
// Calculate Thickness. This is used for the materials
// which thickness is given in % of r.l.
//
/////////////////////////////////////////////////////////
double OraclePixGeoManager::CalculateThickness(double tck,const string& mat) {
  const GeoMaterial* material =  m_pMatMgr->getMaterial(mat);
  double rl = material->getRadLength();
  material->ref();
  material->unref();
  return -1.*rl*tck/100.;
}

int OraclePixGeoManager::moduleType()
{
  int type = 0;
  if (ibl()) {
    if (isBarrel()) {
      type = db()->getInt(m_PixelLayer,"MODULETYPE",m_currentLD);
    }
  } else {
    if(isBarrel()) type = m_currentLD;
    if(isEndcap()) type = m_currentLD+PixelBarrelNLayer();
  }
  return type;
}

int OraclePixGeoManager::moduleType3D()
{
  int type = -1;
  if (!isBarrel()||m_currentLD>0) return type;

  if (ibl()) {
    try {
      type = db()->getInt(m_PixelLayer,"MODULETYPE3D",m_currentLD);
      return type;
    }
    catch(...)
      { 
	return moduleType()+1;
      }
  }
  
  return type;
}



/////////////////////////////////////////////////////////
//
// Si Boards Parameters:
//
/////////////////////////////////////////////////////////
double OraclePixGeoManager::PixelBoardWidth(bool isModule3D) 
{
  if(ibl()&&isModule3D){
    return db()->getDouble(m_PixelModule,"BOARDWIDTH",moduleType3D())*mmcm();
  }

  return db()->getDouble(m_PixelModule,"BOARDWIDTH",moduleType())*mmcm();
}

double OraclePixGeoManager::PixelBoardLength(bool isModule3D) 
{
  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"BOARDLENGTH",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"BOARDLENGTH",moduleType())*mmcm();
}

double OraclePixGeoManager::PixelBoardThickness(bool isModule3D) 
{
  if (m_dc1Geometry && isBarrel() && (m_currentLD == 0)) {
    return 200*Gaudi::Units::micrometer;
  }

  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"BOARDTHICK",moduleType3D())*mmcm();
  return db()->getDouble(m_PixelModule,"BOARDTHICK",moduleType())*mmcm();
}

double  OraclePixGeoManager::PixelBoardActiveLength(bool isModule3D) 
{
  return DesignZActiveArea(isModule3D); 
}


/////////////////////////////////////////////////////////
//
// Hybrid Parameters:
//
/////////////////////////////////////////////////////////
double OraclePixGeoManager::PixelHybridWidth(bool isModule3D) 
{
  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"HYBRIDWIDTH",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"HYBRIDWIDTH",moduleType())*mmcm();
}

double OraclePixGeoManager::PixelHybridLength(bool isModule3D) 
{
  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"HYBRIDLENGTH",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"HYBRIDLENGTH",moduleType())*mmcm();
}

double OraclePixGeoManager::PixelHybridThickness(bool isModule3D) 
{
  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"HYBRIDTHICK",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"HYBRIDTHICK",moduleType())*mmcm();
}
 
/////////////////////////////////////////////////////////
//
// Chip Parameters:
//
/////////////////////////////////////////////////////////

double OraclePixGeoManager::PixelChipWidth(bool isModule3D) 
{
  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"CHIPWIDTH",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"CHIPWIDTH",moduleType())*mmcm();
}

double OraclePixGeoManager::PixelChipLength(bool isModule3D) 
{
  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"CHIPLENGTH",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"CHIPLENGTH",moduleType())*mmcm();
}

double OraclePixGeoManager::PixelChipGap(bool isModule3D) 
{
  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"CHIPGAP",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"CHIPGAP",moduleType())*mmcm();
}

double OraclePixGeoManager::PixelChipOffset(bool isModule3D) 
{
   if(!ibl()||GetLD()!=0||!isBarrel()||!(db()->testField(m_PixelModule,"CHIPOFFSET"))){
     return 0.;
   }

  if(isModule3D)
    return db()->getDouble(m_PixelModule,"CHIPOFFSET",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"CHIPOFFSET",moduleType())*mmcm();
}

double OraclePixGeoManager::PixelChipThickness(bool isModule3D)  {
  if(ibl()&&isModule3D)
    return db()->getDouble(m_PixelModule,"CHIPTHICK",moduleType3D())*mmcm();

  return db()->getDouble(m_PixelModule,"CHIPTHICK",moduleType())*mmcm();
}



/////////////////////////////////////////////////////////
//
// Module services
//
/////////////////////////////////////////////////////////

int OraclePixGeoManager::PixelModuleServiceNumber()
{
  if(!ibl()||GetLD()>0||!isBarrel()) return 0;

  if (db()->getTableSize(m_PixelModuleSvc)) 
    return db()->getTableSize(m_PixelModuleSvc);
  return 0;
}

double OraclePixGeoManager::PixelModuleServiceLength(int svc)
{
  return db()->getDouble(m_PixelModuleSvc,"LENGTH",svc)*mmcm();
}

double OraclePixGeoManager::PixelModuleServiceWidth(int svc)
{
  return db()->getDouble(m_PixelModuleSvc,"WIDTH",svc)*mmcm();
}

double OraclePixGeoManager::PixelModuleServiceThick(int svc)
{
  return db()->getDouble(m_PixelModuleSvc,"THICK",svc)*mmcm();
}

double OraclePixGeoManager::PixelModuleServiceOffsetX(int svc)
{
  return db()->getDouble(m_PixelModuleSvc,"XOFFSET",svc)*mmcm();
}

double OraclePixGeoManager::PixelModuleServiceOffsetY(int svc)
{
  return db()->getDouble(m_PixelModuleSvc,"YOFFSET",svc)*mmcm();
}

double OraclePixGeoManager::PixelModuleServiceOffsetZ(int svc)
{
  return db()->getDouble(m_PixelModuleSvc,"ZOFFSET",svc)*mmcm();
}

int OraclePixGeoManager::PixelModuleServiceFullSize(int svc)
{
  return db()->getInt(m_PixelModuleSvc,"FULLSIZE",svc);
}

int OraclePixGeoManager::PixelModuleServiceModuleType(int svc)
{
  return db()->getInt(m_PixelModuleSvc,"MODULE3D",svc)*mmcm();
}

std::string OraclePixGeoManager::PixelModuleServiceName(int svc)
{
  return db()->getString(m_PixelModuleSvc,"NAME",svc);
}

std::string OraclePixGeoManager::PixelModuleServiceMaterial(int svc)
{
  return db()->getString(m_PixelModuleSvc,"MATERIAL",svc);
}



/////////////////////////////////////////////////////////
//
// Pixel Disks carbon structure
//
/////////////////////////////////////////////////////////
double OraclePixGeoManager::PixelECCarbonRMin(string a) {

  int isup=0;
  if (a == "Inner") {
    isup = 0;
  } else if(a == "Central") {
    isup = 1;
  } else {
    isup = 2;
  }

  return PixelDiskSupportRMin(isup);
}

double OraclePixGeoManager::PixelECCarbonRMax(string a) {
  int isup=0;
  if (a == "Inner") {
    isup = 0;
  } else if(a == "Central") {
    isup = 1;
  } else {
    isup = 2;
  }
  return PixelDiskSupportRMax(isup);
}

double OraclePixGeoManager::PixelECCarbonThickness(string a) {

  int isup=0;
  if (a == "Inner") {
    isup = 0;
  } else if(a == "Central") {
    isup = 1;
  } else {
    isup = 2;
  }
  return PixelDiskSupportThickness(isup);
}


int OraclePixGeoManager::PixelECCarbonMaterialTypeNum(string a) {

  if (dbVersion() < 3) return 0;
  int isup = 0;
  if (a == "Inner") {
    isup = 0;
  } else if(a == "Central") {
    isup = 1;
  } else {
    isup = 2;
  }
  return PixelDiskSupportMaterialTypeNum(isup);
}


/////////////////////////////////////////////////////////
//
// Central Services
//
/////////////////////////////////////////////////////////
//
// This is complicated in the DB...
// If Rmin < 0 this means that the cylinder has to be placed only once
// So I return RMin with its sign, for further processing by the service
// methods.
// If RMax is <0 the thickness is given in % of r.l. which have to be 
// calculated by the method calculatethickness
//
// If Zmin and Zmax have DIFFERENT sign, then the thickness is given in 
// % of r.l....
//

int OraclePixGeoManager::PixelServiceElements(const std::string & type) {
  // FIXME
  /*
  if (dbVersion() < 3) {
    if(isBarrel() ) {
      if(a == "Inside") return (*m_PixelBarrelGeneral)[0]->getInt("NFRAMEIN");
      if(a == "Outside") return (*m_PixelBarrelGeneral)[0]->getInt("NFRAMEOUT");
    }
    if(isEndcap() ) {
      if(a == "Inside") return (*m_PixelEndcapGeneral)[0]->getInt("NFRAMEIN");
      if(a == "Outside") return (*m_PixelEndcapGeneral)[0]->getInt("NFRAMEOUT");   
    }
    return 0;
  } else {
  */
    // a is ignored. Use frame num to distinguish between inside (<1000) and ouside (>=1000). 
  if(type == "simple") return db()->getTableSize(m_PixelSimpleService);
  if(type == "barrel") return db()->getTableSize(m_PixelBarrelService);
  if(type == "endcap") return db()->getTableSize(m_PixelEndcapService);
  if(type == "envelope") return db()->getTableSize(m_PixelEnvelopeService);
  return 0;
  //}
}

// Name used when naming G4/Geo volumes 
std::string OraclePixGeoManager::PixelServiceName(const std::string & type, int index) {

  if (useLegacy() || !getPixelServiceRecordTestField("VOLNAME",type,index)) {
    return "";
  } else {
    return getPixelServiceRecordString("VOLNAME",type,index); 
  }
}


// Flag to say whether volume should be described in both positive and
// negative halves.
bool OraclePixGeoManager::PixelServiceZsymm(const std::string & type, int index) {
  if (dbVersion() < 3 || !getPixelServiceRecordTestField("ZSYMM",type,index)) {
    // If ZSYMM not defined use old logic to determine if volume is
    // duplicated in both positive and negative halves.
    double rmin = getPixelServiceRecordDouble("RIN",type,index);
    double zmin = getPixelServiceRecordDouble("ZIN",type,index);
    return (rmin>0 && zmin > 0.000001); 
  } else {
    return getPixelServiceRecordInt("ZSYMM",type,index);
  }
}


double OraclePixGeoManager::PixelServiceRMin(const std::string & type, int index) {
  return std::abs(getPixelServiceRecordDouble("RIN",type,index)) * mmcm();
}

double OraclePixGeoManager::PixelServiceRMax(const std::string & type, int index) {
  double rtmp =  getPixelServiceRecordDouble("ROUT",type,index);
  // If this is negative this is the thickness of the cyl in % of r.l.
  double rmin = PixelServiceRMin(type,index);
  double rmax = 0;
  if(rtmp > 0) {
    rmax = rtmp * mmcm();
  } else {
    string material = PixelServiceMaterial(type,index);
    rmax = rmin + CalculateThickness(rtmp, material);
  }
  return rmax;
}

double OraclePixGeoManager::PixelServiceRMin2(const std::string & type, int index) {
  if (!getPixelServiceRecordTestField("RIN2",type,index)) {
    return 0;
  } else {
    return getPixelServiceRecordDouble("RIN2",type,index) * Gaudi::Units::mm;
  }
}

double OraclePixGeoManager::PixelServiceRMax2(const std::string & type, int index) {
  if (!getPixelServiceRecordTestField("ROUT2",type,index)) {
    return 0;
  } else {
    return getPixelServiceRecordDouble("ROUT2",type,index) * Gaudi::Units::mm;
  }
}

double OraclePixGeoManager::PixelServiceZMin(const std::string & type, int index) {
  return getPixelServiceRecordDouble("ZIN",type,index) * mmcm();
}

double OraclePixGeoManager::PixelServiceZMax(const std::string & type, int index) {
  double zout =  getPixelServiceRecordDouble("ZOUT",type,index);
  double zmin  =  PixelServiceZMin(type,index);
  double zmax = 0;
  // If zin and zout are opposite sign then zout is the thickness of the cyl in % of r.l.
  if(zmin * zout > -0.000000001) {
    zmax = zout * mmcm();
  } else {
    string material = PixelServiceMaterial(type,index);
    double sign = (zmin > 0) ? 1: -1;
    zmax = zmin + sign*CalculateThickness(zout, material);
  }
  return zmax;
}

double OraclePixGeoManager::PixelServicePhiLoc(const std::string & type, int index) {
  if (!getPixelServiceRecordTestField("PHI",type,index)) {
    return 0;
  } else {
    return getPixelServiceRecordDouble("PHI",type,index) * Gaudi::Units::degree; 
  }
}

double OraclePixGeoManager::PixelServiceWidth(const std::string & type, int index) {
  if (!getPixelServiceRecordTestField("WIDTH",type,index)) {
    return 0;
  } else {
    // Can be in degree or Gaudi::Units::mm. Leave it to GeoPixelServices to interpret.    
    return getPixelServiceRecordDouble("WIDTH",type,index); 
  }
}

int OraclePixGeoManager::PixelServiceRepeat(const std::string & type, int index) {
  if (!getPixelServiceRecordTestField("REPEAT",type,index)) {
    return 0;
  } else {
    return getPixelServiceRecordInt("REPEAT",type,index); 
  }
}

std::string OraclePixGeoManager::PixelServiceShape(const std::string & type, int index) {
  if (type == "simple") return "TUBE";
  if (!getPixelServiceRecordTestField("SHAPE",type,index)) {
    return "TUBE";
  } else {
    return getPixelServiceRecordString("SHAPE",type,index); 
  }
}


int OraclePixGeoManager::PixelServiceShift(const std::string & type, int index) {
  if (!getPixelServiceRecordTestField("SHIFT",type,index)) {
    return 0;
  } else {
    return getPixelServiceRecordInt("SHIFT",type,index); 
  }
}


int OraclePixGeoManager::PixelServiceLD(const std::string & type, int index) {
  return getPixelServiceRecordInt("LAYERNUM",type,index)-1;
}

string OraclePixGeoManager::PixelServiceMaterial(const std::string & type, int index) {

  int imat = 0;
  if (type != "simple") {
    imat = getPixelServiceRecordInt("MATERIAL",type,index);
  }
  std::string materialName;
  if (!imat) {
    materialName = getPixelServiceRecordString("MATERIALNAME",type,index);
  } else {
    // Old
    if(type == "barrel") {
      string mat[11] = {
	"std::Berillia",
	"std::Carbon",
	"pix::Services",
	"pix::Titanium",
	"pix::MatPP11",
	"pix::MatPP12",
	"pix::MatPP13",
	"pix::MatPP14",
	"pix::MatPP15",
	"pix::MatPP16",
	"pix::MatPP17"};
      materialName =  mat[imat-1];
    } 
    if(type == "endcap") {
      string mat[4] = {
	"std::Berillia",
	"std::Carbon",
	"pix::ECServices",
	"pix::Disk"};
      materialName =  mat[imat-1];
    }
  }
  return materialName;
}


int OraclePixGeoManager::PixelServiceFrameNum(const std::string & type, int index) {
  // In older version frame num indicated "inside" or "outside"
  // 0-999:  Inside
  // >=1000: Outside
  // No frame number in simple table.
  if (type == "simple") return index+1;
  int framenum = getPixelServiceRecordInt("FRAMENUM",type,index);
  if (framenum <= 0) return index+1;
  if (dbVersion() < 3) {
    if(type == "barrel") {
      if (index >= m_barrelInFrames) framenum += 1000;
    }
    if(type == "endcap") {
      if (index >= m_endcapInFrames) framenum += 1000;
    }
  }
  return framenum;
  // FIXME
  /*
  if (dbVersion() < 3) return framenum;
  if (type == "Inside" && framenum < 1000) return framenum;
  if (type == "Outside" && framenum >= 1000) return framenum%1000;
  return -1; 
  */
}

// Access child/envelope service parameters
int OraclePixGeoManager::PixelServiceEnvelopeNum(const std::string & type, int index) {

  if (type != "envelope") return 0;

  try{
    int envnum = getPixelServiceRecordInt("ENVNUM",type,index);
    return envnum;
  }
  catch(...)
    {}

  return 0;
}

int OraclePixGeoManager::PixelServiceParentEnvelopeNum(const std::string & type, int index) {

  //  if (type == "envelope") return 0;

  if (type == "envelope"){
    try{
      int envnum = getPixelServiceRecordInt("ENVPARENT",type,index);
      return envnum;
    }
    catch(...)
      {}
  }
  else {

    try{
      int envnum = getPixelServiceRecordInt("ENVELOPE",type,index);
      return envnum;
    }
    catch(...)
      {
      }
  }


  return 0;
}

std::string OraclePixGeoManager::getPixelServiceRecordString(const std::string & name, const std::string & type, int index) {
  IRDBRecordset_ptr recordSet = getPixelServiceRecordset(type);
  return db()->getString(recordSet, name, index);
}

int OraclePixGeoManager::getPixelServiceRecordInt(const std::string & name, const std::string & type, int index) {
  IRDBRecordset_ptr recordSet = getPixelServiceRecordset(type);
  return db()->getInt(recordSet, name, index);
}


double OraclePixGeoManager::getPixelServiceRecordDouble(const std::string & name, const std::string & type, int index) {
  IRDBRecordset_ptr recordSet = getPixelServiceRecordset(type);
  return db()->getDouble(recordSet, name, index);
}

bool OraclePixGeoManager::getPixelServiceRecordTestField(const std::string & name, const std::string & type, int index) {
  try {
    IRDBRecordset_ptr recordSet = getPixelServiceRecordset(type);
    return db()->testField(recordSet, name, index);
  }
  catch(...){}
  return false;
}


// Returns IRDBRecordset
IRDBRecordset_ptr  OraclePixGeoManager::getPixelServiceRecordset(const std::string & type) {
  // m_barrelInFrames and m_endcapInFrames should be zero in dbVersion >= 3
  IRDBRecordset_ptr recordSet;
  if (type == "simple") {
    recordSet = m_PixelSimpleService;
  } else if(type == "barrel") {
    recordSet = m_PixelBarrelService;
    //if(type != "Inside") index += m_barrelInFrames;
  } else if(type == "endcap") {
    recordSet = m_PixelEndcapService;
    //if(type != "Inside") index += m_endcapInFrames;
  } else if(type == "envelope") {
    recordSet = m_PixelEnvelopeService;
    //if(type != "Inside") index += m_endcapInFrames;
  } else {
    msg(MSG::ERROR) << "ERROR:  getPixelServiceRecord(), neither Barrel of Endcap selected!" << endmsg;
  }
  return recordSet;
} 

double OraclePixGeoManager::PixelECCablesThickness() 
{
  double tck = db()->getDouble(m_PixelDisk,"CABLETHICK",m_currentLD);
  if( tck > 0.) {
    return tck*mmcm();
  } else {    
    std::string matName =  getMaterialName("DiskCable", m_currentLD);
    return CalculateThickness(tck,matName);
  }
}

int 
OraclePixGeoManager::PixelCableElements()
{
  if (dbVersion() < 3) return m_legacyManager->PixelCableElements();
  return db()->getTableSize(m_PixelBarrelCable);
}

int 
OraclePixGeoManager::PixelCableLayerNum(int index)
{
  if (dbVersion() < 3) return 0;
  return db()->getInt(m_PixelBarrelCable,"LAYER",index);
}

int 
OraclePixGeoManager::PixelCableBiStaveNum(int index)
{
  if (dbVersion() < 3) return 0;
  return db()->getInt(m_PixelBarrelCable,"BISTAVE",index);
}


double 
OraclePixGeoManager::PixelCableZStart(int index)
{
  if (dbVersion() < 3) return m_legacyManager->PixelCableZStart(index);
  return db()->getDouble(m_PixelBarrelCable,"ZSTART",index) * Gaudi::Units::mm;
}

double 
OraclePixGeoManager::PixelCableZEnd(int index)
{
  if (dbVersion() < 3) return m_legacyManager->PixelCableZEnd(index);
  return db()->getDouble(m_PixelBarrelCable,"ZEND",index) * Gaudi::Units::mm;
}

double 
OraclePixGeoManager::PixelCableWidth(int index)
{
  if (dbVersion() < 3) return m_legacyManager->PixelCableWidth(index);
  return db()->getDouble(m_PixelBarrelCable,"WIDTH",index) * Gaudi::Units::mm;
}

double 
OraclePixGeoManager::PixelCableThickness(int index)
{
  if (dbVersion() < 3) return m_legacyManager->PixelCableThickness(index);
  return db()->getDouble(m_PixelBarrelCable,"THICK",index) * Gaudi::Units::mm;
}

double 
OraclePixGeoManager::PixelCableStackOffset(int index)
{
  if (dbVersion() < 3) return m_legacyManager->PixelCableStackOffset(index);
  return db()->getDouble(m_PixelBarrelCable,"STACKPOS",index) * Gaudi::Units::mm;
}

double 
OraclePixGeoManager::PixelCableWeight(int index)
{
  if (dbVersion() < 3) return 0;
  return db()->getDouble(m_PixelBarrelCable,"WEIGHT",index) * GeoModelKernelUnits::g;
}

std::string
OraclePixGeoManager::PixelCableLabel(int index)
{
  if (dbVersion() < 3) return m_legacyManager->PixelCableLabel(index);
  return db()->getString(m_PixelBarrelCable,"LABEL",index);
}


//
// Version of the Geometry
//

int OraclePixGeoManager::determineDbVersion() {
  // This determines a version depending on various changes in the database;
  int version = 0;

  if (!(*m_PixelLayer)[0]->isFieldNull("PHIOFMODULEZERO")) version = 1;
  if (m_PixelReadout->size() != 0) version = 2;
  if (m_weightTable->size() != 0) version = 3;
  if (m_PixelTMT->size() != 0) version = 4; // Removed all legacy tables

  return version;
}



std::string OraclePixGeoManager::getMaterialName(const std::string & volumeName, int layerdisk, int typenum) {
  return m_materialMap->getMaterial(layerdisk, typenum, volumeName);
}


void OraclePixGeoManager::addDefaultMaterials() {
  // This is for backward compatibilty. Newer geometies get the
  // gets them from the database.
  m_materialMap->addMaterial(0,0,"Sensor","std::Silicon");
  m_materialMap->addMaterial(0,0,"Chip","pix::Chip");
  m_materialMap->addMaterial(0,0,"Hybrid","pix::Hybrid");
  m_materialMap->addMaterial(0,0,"Omega","pix::MatOmega");
  m_materialMap->addMaterial(0,0,"AlTube","pix::MatAlTube");
  m_materialMap->addMaterial(1,0,"AlTube","pix::MatAlTube");
  m_materialMap->addMaterial(2,0,"AlTube","pix::MatAlTubeFix");
  m_materialMap->addMaterial(0,1,"Fluid","pix::MatCap1");
  m_materialMap->addMaterial(0,2,"Fluid","pix::MatCap2");
  m_materialMap->addMaterial(0,0,"TMT","pix::MatTMT");
  m_materialMap->addMaterial(0,0,"GlueOmegaStave","pix::MatGlue");
  m_materialMap->addMaterial(0,0,"Connector","pix::MatConn");
  m_materialMap->addMaterial(0,0,"PigtailCyl","pix::MatPigtail");
  m_materialMap->addMaterial(0,0,"PigtailFlat","pix::MatPigtail");
  m_materialMap->addMaterial(0,0,"Cable","pix::MatT0");
  m_materialMap->addMaterial(0,0,"DiskCable","pix::ECCables");
  m_materialMap->addMaterial(0,0,"DiskSupport","pix::Disk");
  m_materialMap->addMaterial(0,0,"Frame","std::Carbon");
  m_materialMap->addMaterial(0,0,"EndCone","std::Carbon");
}

std::string OraclePixGeoManager::getLD_Label()
{
  std::ostringstream o;
  if(isBarrel()) {
     if (m_currentLD == 0) {
       o << "BL";
     } else {
       o << "L" << m_currentLD;
     }
  } else {
    o << "D" << m_currentLD;
  }
  return o.str();
}
 
int OraclePixGeoManager::PixelBarrelMajorVersion()
{ 
  return static_cast<int>((*m_PixelBarrelGeneral)[0]->getDouble("VERSION"));
}

int OraclePixGeoManager::PixelBarrelMinorVersion()
{ 
  return static_cast<int>(((*m_PixelBarrelGeneral)[0]->getDouble("VERSION") - PixelBarrelMajorVersion())*10 + 0.5);
}

int OraclePixGeoManager::PixelEndcapMajorVersion()
{
  return static_cast<int>((*m_PixelEndcapGeneral)[0]->getDouble("VERSION"));
}

int OraclePixGeoManager::PixelEndcapMinorVersion() 
{ 
  return static_cast<int>(((*m_PixelEndcapGeneral)[0]->getDouble("VERSION") - PixelEndcapMajorVersion())*10 + 0.5);
}


std::string OraclePixGeoManager::versionDescription() const
{
  std::string description;
  if (db()->testField(m_PixelSwitches,"DESCRIPTION")) {
    description = db()->getString(m_PixelSwitches,"DESCRIPTION");
  }
  return description;
}

std::string OraclePixGeoManager::versionName() const
{
  std::string name;
  if (db()->testField(m_PixelSwitches,"VERSIONNAME")) {
    name = db()->getString(m_PixelSwitches,"VERSIONNAME");
  }
  return name;
}

std::string OraclePixGeoManager::versionLayout() const
{
  std::string layout;
  if (db()->testField(m_PixelSwitches,"LAYOUT")) {
    layout = db()->getString(m_PixelSwitches,"LAYOUT");
  }
  return layout;
}


double OraclePixGeoManager::PixelRMin() 
{
  if (db()->getTableSize(m_PixelEnvelope)) {
    double rmin = PixelEnvelopeRMin(0);  
    for (unsigned int i = 1; i < db()->getTableSize(m_PixelEnvelope); i++) {
      rmin = std::min(rmin, PixelEnvelopeRMin(i));
    } 
    return rmin;
  } else {      
    return db()->getDouble(m_PixelCommon,"RMIN")*mmcm();
  }
}

double OraclePixGeoManager::PixelRMax() 
{
  if (db()->getTableSize(m_PixelEnvelope)) {
    double  rmax = PixelEnvelopeRMax(0);  
    for (unsigned int i = 1; i < db()->getTableSize(m_PixelEnvelope); i++) {
      rmax = std::max(rmax, PixelEnvelopeRMax(i));
    } 
    return rmax;
  } else {      
    return db()->getDouble(m_PixelCommon,"RMAX")*mmcm();
  }
}

double OraclePixGeoManager::PixelHalfLength() 
{

  if (db()->getTableSize(m_PixelEnvelope)) {
    // The table should contain only +ve z values.
    return PixelEnvelopeZ(db()->getTableSize(m_PixelEnvelope) - 1);
  } else {
    return db()->getDouble(m_PixelCommon,"HALFLENGTH")*mmcm();
  }
}

bool OraclePixGeoManager::PixelSimpleEnvelope() 
{
  // Return true if the envelope can be built as a simple tube.
  // otherwise it will be built as a PCON.
  // True if size is 0 or 1.
  return (!(db()->getTableSize(m_PixelEnvelope) > 1));
}

unsigned int OraclePixGeoManager::PixelEnvelopeNumPlanes() 
{
  return db()->getTableSize(m_PixelEnvelope);
}

double OraclePixGeoManager::PixelEnvelopeZ(int i) 
{
  double zmin =  db()->getDouble(m_PixelEnvelope,"Z",i) * Gaudi::Units::mm;
  if (zmin < 0) msg(MSG::ERROR) << "PixelEnvelope table should only contain +ve z values" << endmsg;
  return std::abs(zmin);
}

double OraclePixGeoManager::PixelEnvelopeRMin(int i) 
{
  return db()->getDouble(m_PixelEnvelope,"RMIN",i) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelEnvelopeRMax(int i) 
{
  return db()->getDouble(m_PixelEnvelope,"RMAX",i) * Gaudi::Units::mm;
}


int OraclePixGeoManager::PixelBarrelNLayer() 
{
  return db()->getInt(m_PixelBarrelGeneral,"NLAYER");
}

// m_PixelBarrelGeneral
double OraclePixGeoManager::PixelBarrelRMin() 
{
  return db()->getDouble(m_PixelBarrelGeneral,"RMIN")*mmcm();
}

double OraclePixGeoManager::PixelBarrelRMax() 
{
  return db()->getDouble(m_PixelBarrelGeneral,"RMAX")*mmcm();
}

double OraclePixGeoManager::PixelBarrelHalfLength() 
{
  return db()->getDouble(m_PixelBarrelGeneral,"HALFLENGTH")*mmcm();
}

// Described in general services for later geometries.
bool OraclePixGeoManager::oldFrame()
{
  if (useLegacy()) return m_legacyManager->oldFrame();
  return false;
} 

// For new geometry a detailed frame is built.
bool OraclePixGeoManager::detailedFrame()
{
  return db()->getTableSize(m_PixelFrame);
}
  
int OraclePixGeoManager::PixelFrameSections()
{
  return db()->getTableSize(m_PixelFrame);
}

double OraclePixGeoManager::PixelFrameRMinSide(int sectionIndex)
{
  return db()->getDouble(m_PixelFrame, "RMINSIDE", sectionIndex) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFrameRMaxSide(int sectionIndex)
{
  return db()->getDouble(m_PixelFrame, "RMAXSIDE", sectionIndex) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFrameSideWidth(int sectionIndex)
{
  return db()->getDouble(m_PixelFrame, "SIDEWIDTH", sectionIndex) * Gaudi::Units::mm;
} 
 
double OraclePixGeoManager::PixelFrameZMin(int sectionIndex)
{ 
  return db()->getDouble(m_PixelFrame, "ZMIN", sectionIndex) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFrameZMax(int sectionIndex)
{ 
  return db()->getDouble(m_PixelFrame, "ZMAX", sectionIndex) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFramePhiStart(int sectionIndex)
{
  return db()->getDouble(m_PixelFrame, "PHISTART", sectionIndex) * Gaudi::Units::deg;
}
 
int OraclePixGeoManager::PixelFrameNumSides(int sectionIndex)
{
  return db()->getInt(m_PixelFrame, "NUMSIDES", sectionIndex);
}

bool OraclePixGeoManager::PixelFrameMirrorSides(int sectionIndex)
{
  return db()->getInt(m_PixelFrame, "MIRRORSIDES", sectionIndex);
}
									
std::string OraclePixGeoManager::PixelFrameSideMaterial(int sectionIndex)
{
  return db()->getString(m_PixelFrame, "SIDEMATERIAL", sectionIndex);
}

std::string OraclePixGeoManager::PixelFrameCornerMaterial(int sectionIndex)
{
  return db()->getString(m_PixelFrame, "CORNERMATERIAL", sectionIndex);
} 

int OraclePixGeoManager::PixelFrameSectionFromIndex(int sectionIndex)
{
  return db()->getInt(m_PixelFrame,"SECTION",sectionIndex);
}
  
void 
OraclePixGeoManager::makeFrameIndexMap()
{
  if (!m_frameElementMap) {
    m_frameElementMap = new std::map<int,std::vector<int> >;
    for (unsigned int i = 0; i < db()->getTableSize(m_PixelFrameSect); ++i) {
      int section = db()->getInt(m_PixelFrameSect,"SECTION",i);
      (*m_frameElementMap)[section].push_back(i);
    }
  }
}

int OraclePixGeoManager::getFrameElementIndex(int sectionIndex, int element)
{
  // make map if it is not already made.
  makeFrameIndexMap();

  int section = PixelFrameSectionFromIndex(sectionIndex);

  int newIndex = -1;
  std::map<int,std::vector<int> >::const_iterator iter = m_frameElementMap->find(section);
  if (iter ==  m_frameElementMap->end()) {
    // Should never be the case as PixelFrameNumSideElements should generally be called first
    msg(MSG::ERROR) << "Frame section " << section << " has no elements." << endmsg;
  } else {
    const std::vector<int> & vec = iter->second;
    if (static_cast<unsigned int>(element) >= vec.size()) {
      msg(MSG::ERROR) << "Element index " << element << " for section " << section << " out of range." << endmsg;
    } else {
      newIndex = vec[element];
    }
  }
  return newIndex;
}


int OraclePixGeoManager::PixelFrameNumSideElements(int sectionIndex)
{ 
  // make map if it is not already made.
  makeFrameIndexMap();

  int section = PixelFrameSectionFromIndex(sectionIndex);
  int numElements = 0;

  std::map<int,std::vector<int> >::const_iterator iter = m_frameElementMap->find(section);
  if (iter ==  m_frameElementMap->end()) {
    msg(MSG::DEBUG) << "Frame section " << section << " has no elements." << endmsg;
  } else {
    numElements = iter->second.size();
  }
  return numElements;
}
 
double OraclePixGeoManager::PixelFrameElementZMin1(int sectionIndex, int element)
{
  int index = getFrameElementIndex(sectionIndex, element);
  if (index < 0) return 0; // Error message already printed in getFrameElementIndex.
  return db()->getDouble(m_PixelFrameSect, "ZMIN1", index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFrameElementZMin2(int sectionIndex, int element)
{
  int index = getFrameElementIndex(sectionIndex, element);
  if (index < 0) return 0; // Error message already printed in getFrameElementIndex.
  return db()->getDouble(m_PixelFrameSect, "ZMIN2", index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFrameElementZMax1(int sectionIndex, int element)
{
  int index = getFrameElementIndex(sectionIndex, element);
  if (index < 0) return 0; // Error message already printed in getFrameElementIndex.
  return db()->getDouble(m_PixelFrameSect, "ZMAX1", index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFrameElementZMax2(int sectionIndex, int element)
{
  int index = getFrameElementIndex(sectionIndex, element);
  if (index < 0) return 0; // Error message already printed in getFrameElementIndex.
  return db()->getDouble(m_PixelFrameSect, "ZMAX2", index) * Gaudi::Units::mm;
}

int OraclePixGeoManager::PixelStaveIndex(int layer)
{
  if (!ibl()) return 0;
  if (!db()->testField(m_PixelLayer,"STAVEINDEX",layer)) return 0;
  return db()->getInt(m_PixelLayer,"STAVEINDEX",layer);
}

int OraclePixGeoManager::PixelStaveLayout()
{
  if (!ibl()) return 0;
  int defaultLayout = 0;
  int index = PixelStaveIndex(m_currentLD);

  //  if(m_currentLD==0)
  //    std::cout<<"OraclePixGeoManager::PixelStaveLayout() : "<<m_currentLD<<" "<<index<<" "<<db()->getInt(m_PixelStave,"LAYOUT",index)<<std::endl;


  if (!db()->testField(m_PixelStave,"LAYOUT",index)) return defaultLayout;
  return db()->getInt(m_PixelStave,"LAYOUT",index);
}

int OraclePixGeoManager::PixelStaveAxe()
{
  if (!ibl()) return 0;
  int index = PixelStaveIndex(m_currentLD);

  if (db()->testField(m_PixelStave,"STAVEAXE",index))
    return db()->getInt(m_PixelStave,"STAVEAXE",index);
  return 0;
}

double OraclePixGeoManager::PixelLayerRadius() 
{
  double radius = db()->getDouble(m_PixelLayer,"RLAYER",m_currentLD)*mmcm();
  if (msgLvl(MSG::DEBUG)) msg(MSG::DEBUG) << "PixelLayerRadius for layer " << m_currentLD
      << " is " << radius
      << endmsg;
  return radius;
}

double OraclePixGeoManager::PixelLayerGlobalShift() 
{
  if (db()->testField(m_PixelLayer,"GBLSHIFT",m_currentLD))
    return db()->getDouble(m_PixelLayer,"GBLSHIFT",m_currentLD);
  return 0.;
}

double OraclePixGeoManager::PixelLadderLength() 
{
  if (useLegacy()) return m_legacyManager->PixelLadderLength(); 
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"ENVLENGTH",index)*Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelLadderWidthClearance() 
{
  if (useLegacy()) return 0.9*Gaudi::Units::mm; 
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"CLEARANCEY",index)*Gaudi::Units::mm;  
}

// Only used if ladder thickness is automatically calculated it, ie ENVTHICK = 0
// IBL only
double OraclePixGeoManager::PixelLadderThicknessClearance() 
{
  int index = PixelStaveIndex(m_currentLD);
  if (db()->testField(m_PixelStave,"CLEARANCEX",index)) {
    return db()->getDouble(m_PixelStave,"CLEARANCEX",index)*Gaudi::Units::mm;  
  }
  return 0.1*Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelLadderThickness() 
{
  if (useLegacy()) return m_legacyManager->PixelLadderThickness();  // 2*1.48972 mm
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"ENVTHICK",index)*Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelLadderTilt() 
{
  return db()->getDouble(m_PixelLayer,"STAVETILT",m_currentLD)*Gaudi::Units::deg;
}

double OraclePixGeoManager::PixelLadderServicesX() 
{
  if (useLegacy()) return m_legacyManager->PixelLadderServicesX(); // 1.48972 mm
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"SERVICEOFFSETX",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelLadderServicesY() 
{
  if (useLegacy()) return m_legacyManager->PixelLadderServicesY();  // 3mm
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"SERVICEOFFSETY",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelLadderCableOffsetX() 
{
  if (useLegacy()) return m_legacyManager->PixelLadderCableOffsetX(); // 0
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"CABLEOFFSETX",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelLadderCableOffsetY() 
{
  if (useLegacy()) return m_legacyManager->PixelLadderCableOffsetY();  // 4mm
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"CABLEOFFSETY",index) * Gaudi::Units::mm;
}

// IBL only
double OraclePixGeoManager::PixelLadderSupportThickness() 
{
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"SUPPORTTHICK",index) * Gaudi::Units::mm;
}

// IBL only
double OraclePixGeoManager::PixelLadderSupportWidth() 
{
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"SUPPORTWIDTH",index) * Gaudi::Units::mm;
}





// IBL only
double OraclePixGeoManager::PixelLadderBentStaveAngle() 
{
  if (!db()->testFieldTxt(m_PixelConicalStave,"BENTSTAVEANGLE")) return 0;
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelConicalStave,"BENTSTAVEANGLE",index);
}

// IBL only
int OraclePixGeoManager::PixelBentStaveNModule() 
{
  if (!db()->testFieldTxt(m_PixelConicalStave,"BENTSTAVENMODULE")) return 0;
  int index = PixelStaveIndex(m_currentLD);
  return db()->getInt(m_PixelConicalStave,"BENTSTAVENMODULE",index);
}

double OraclePixGeoManager::PixelLadderModuleDeltaZ()
{
  int index = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"MODULEDZ",index);
}

// IBL only
double OraclePixGeoManager::PixelLadderSupportLength() 
{
  int index = PixelStaveIndex(m_currentLD);
  if (db()->testField(m_PixelStave,"SUPPORTHLENGTH",index)) {
    double halflength = db()->getDouble(m_PixelStave,"SUPPORTHLENGTH",index) * Gaudi::Units::mm;
    if (halflength > 0)  return 2 * halflength;
  } 
  double safety = 0.01*Gaudi::Units::mm;
  return PixelLadderLength() - safety;
}

// IBL detailed stave support only

GeoTrf::Vector3D OraclePixGeoManager::IBLStaveRotationAxis() 
{
  // set layer to 0  (in order to read read IBL data)
  int currentLD_tmp = m_currentLD;
  m_currentLD = 0;

  double boardThick = PixelBoardThickness();
  double chipThick = PixelChipThickness();
  double chipGap = PixelChipGap();

  double xCenterCoolingPipe = boardThick*.5+chipThick+chipGap+                // from sensor sensor to plate
    IBLStaveFacePlateThickness() + IBLStaveFacePlateGreaseThickness() +       // plate thickness (grease + plate)
    IBLStaveTubeMiddlePos();                                                  // from plate to colling pipe center
  double yCenterCoolingPipe = IBLStaveMechanicalStaveOffset();
  GeoTrf::Vector3D centerCoolingPipe(xCenterCoolingPipe, yCenterCoolingPipe, 0.);

  m_currentLD = currentLD_tmp;  
  return centerCoolingPipe;
}


double OraclePixGeoManager::IBLStaveRadius() 
{
  // set layer to 0  (in order to read read IBL data)
  int currentLD_tmp = m_currentLD;
  m_currentLD = 0;

  //  Point that defines the center of the cooling pipe
  GeoTrf::Vector3D centerCoolingPipe_inv = -IBLStaveRotationAxis();
  GeoTrf::Vector3D origin(0.,0.,0.);
  double layerRadius = PixelLayerRadius();
  double ladderTilt  = PixelLadderTilt();
  
  // Transforms
  GeoTrf::Transform3D staveTrf = GeoTrf::RotateZ3D(ladderTilt)*GeoTrf::Translate3D(centerCoolingPipe_inv.x(),centerCoolingPipe_inv.y(),centerCoolingPipe_inv.z());
  GeoTrf::Vector3D sensorPos = staveTrf*origin;
  
  double yPos = sensorPos.y();
  GeoTrf::Vector3D sensorPos_layer(sqrt(layerRadius*layerRadius-yPos*yPos),yPos,0.);
  
  double staveRadius = sensorPos_layer.x()-sensorPos.x();

  m_currentLD = currentLD_tmp;  
  return staveRadius;
}


double OraclePixGeoManager::IBLStaveFacePlateThickness() 
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"FACEPLATETHICK",index)) {
    double thickness = db()->getDouble(m_PixelIBLStave,"FACEPLATETHICK",index) * Gaudi::Units::mm;
    if (thickness > 0)  return thickness ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveMechanicalStaveWidth()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"STAVEWIDTH",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"STAVEWIDTH",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveMechanicalStaveEndBlockLength()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"ENDBLOCKLENGTH",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"ENDBLOCKLENGTH",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveMechanicalStaveEndBlockFixPoint()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"ENDBLOCKFIXINGPOS",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"ENDBLOCKFIXINGPOS",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveMechanicalStaveEndBlockOmegaOverlap()
{
  // try and catch (param availbale only if db tag > IBL-03-00-00)
  try{
    int index=0;
    if (db()->testField(m_PixelIBLStave,"ENDBLOCKOMEGAOVERLAP",index)) {
      double value = db()->getDouble(m_PixelIBLStave,"ENDBLOCKOMEGAOVERLAP",index) * Gaudi::Units::mm;
      return value ;
    } 
    return 0.0;
  }
  catch(...){}
  return 0.;
}

double OraclePixGeoManager::IBLStaveLength()
{
 // try and catch (param availbale only if db tag > IBL-03-00-00)
  try
    {
      int index=0;
      if (db()->testField(m_PixelIBLStave,"STAVELENGTH",index)) {
	double value = db()->getDouble(m_PixelIBLStave,"STAVELENGTH",index) * Gaudi::Units::mm;
	return value ;
      } 
    }
  catch(...)
    {
      // FIXME : patch for initial IBL geometry (SES)
      //           IBL stave length not eqal to other stave length 
    }  
  
  return 748.0 * Gaudi::Units::mm;  
}

double OraclePixGeoManager:: IBLStaveMechanicalStaveOffset(bool isModule3D)
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (!isModule3D&&db()->testField(m_PixelIBLStave,"MODULELATERALOFFSET",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"MODULELATERALOFFSET",index) * Gaudi::Units::mm;
    return value ;
  } 
  if (isModule3D&&db()->testField(m_PixelIBLStave,"MODULELATERALOFFSET3D",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"MODULELATERALOFFSET3D",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveMechanicalStaveModuleOffset()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"STAVETOMODULEGAP",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"STAVETOMODULEGAP",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveTubeOuterDiameter()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"TUBEOUTERDIAM",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"TUBEOUTERDIAM",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveTubeInnerDiameter()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"TUBEINNERDIAM",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"TUBEINNERDIAM",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveTubeMiddlePos()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"TUBEMIDDLEPOS",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"TUBEMIDDLEPOS",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveFlexLayerThickness()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"FLEXLAYERTHICK",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"FLEXLAYERTHICK",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveFlexBaseThickness()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"FLEXBASETHICK",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"FLEXBASETHICK",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveFlexWidth()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"FLEXWIDTH",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"FLEXWIDTH",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLStaveFlexOffset()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"FLEXOFFSET",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"FLEXOFFSET",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}


double OraclePixGeoManager::IBLStaveOmegaThickness()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"OMEGATHICK",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"OMEGATHICK",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLStaveOmegaEndCenterX()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"OMEGAENDCENTERX",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"OMEGAENDCENTERX",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}
double OraclePixGeoManager::IBLStaveOmegaEndCenterY()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"OMEGAENDCENTERY",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"OMEGAENDCENTERY",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}
double OraclePixGeoManager::IBLStaveOmegaEndRadius()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"OMEGAENDRADIUS",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"OMEGAENDRADIUS",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}
double OraclePixGeoManager::IBLStaveOmegaEndAngle()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"OMEGAENDANGLE",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"OMEGAENDANGLE",index) * Gaudi::Units::deg;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLStaveOmegaMidCenterX()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"OMEGAMIDCENTERX",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"OMEGAMIDCENTERX",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLStaveOmegaMidRadius()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"OMEGAMIDRADIUS",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"OMEGAMIDRADIUS",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}
double OraclePixGeoManager::IBLStaveOmegaMidAngle()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"OMEGAOPENINGANGLE",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"OMEGAOPENINGANGLE",index) * Gaudi::Units::deg;
    return value ;
  } 
  return 0.0;
}

int OraclePixGeoManager::IBLStaveModuleNumber_AllPlanar()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"NMODULE",index)) {
    int value = db()->getInt(m_PixelIBLStave,"NMODULE",index);
    if (value > 0)  return value ;
  } 
  return 0;
}

int OraclePixGeoManager::IBLStaveModuleNumber()
{
  return m_PlanarModuleNumber+m_3DModuleNumber;

}

double OraclePixGeoManager::IBLStaveModuleGap()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"MODULEGAP",index)) {
    double value = db()->getDouble(m_PixelIBLStave,"MODULEGAP",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

int OraclePixGeoManager::IBLStaveModuleType() 
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLStave,"MODULETYPE",index)) {
    int value = db()->getInt(m_PixelIBLStave,"MODULETYPE",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0;
}

double OraclePixGeoManager::IBLStaveFacePlateGreaseThickness()
{
  // try and catch (param availbale only if db tag > IBL-03-00-00)
  try{
    int index=0;
    if (db()->testField(m_PixelIBLGlueGrease,"FACEPLATEGREASETHICK",index)) {
      double value = db()->getDouble(m_PixelIBLGlueGrease,"FACEPLATEGREASETHICK",index) * Gaudi::Units::mm;
      return value ;
    }
    return 0.;
  }
  catch(...){}
  return 0.;
}

double OraclePixGeoManager::IBLStaveFacePlateGlueThickness()
{
  // try and catch (param availbale only if db tag > IBL-03-00-00)
  try{
    int index=0;
    if (db()->testField(m_PixelIBLGlueGrease,"FACEPLATEGLUETHICK",index)) {
      double value = db()->getDouble(m_PixelIBLGlueGrease,"FACEPLATEGLUETHICK",index) * Gaudi::Units::mm;
      return value ;
    }
    return 0.;
  }
  catch(...) {}
  return 0.;
}

double OraclePixGeoManager::IBLStaveTubeGlueThickness()
{
  // try and catch (param availbale only if db tag > IBL-03-00-00)
  try{
    int index=0;
    if (db()->testField(m_PixelIBLGlueGrease,"TUBEGLUETHICK",index)) {
      double value = db()->getDouble(m_PixelIBLGlueGrease,"TUBEGLUETHICK",index) * Gaudi::Units::mm;
      return value ;
    }
    return 0.;
  }
  catch(...) {}
  return 0.;
}

double OraclePixGeoManager::IBLStaveOmegaGlueThickness()
{
  // try and catch (param availbale only if db tag > IBL-03-00-00)
  try{
    int index=0;
    if (db()->testField(m_PixelIBLGlueGrease,"OMEGAGLUETHICK",index)) {
      double value = db()->getDouble(m_PixelIBLGlueGrease,"OMEGAGLUETHICK",index) * Gaudi::Units::mm;
      return value ;
    }
    return 0.;
  }
  catch(...){}
  return 0.;
}


double OraclePixGeoManager:: IBLSupportRingWidth()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLSupport,"STAVERINGWIDTH",index)) {
    double value = db()->getDouble(m_PixelIBLSupport,"STAVERINGWIDTH",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLSupportRingInnerRadius()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLSupport,"STAVERINGINNERRADIUS",index)) {
    double value = db()->getDouble(m_PixelIBLSupport,"STAVERINGINNERRADIUS",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLSupportRingOuterRadius()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLSupport,"STAVERINGOUTERRADIUS",index)) {
    double value = db()->getDouble(m_PixelIBLSupport,"STAVERINGOUTERRADIUS",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}


double OraclePixGeoManager:: IBLSupportMechanicalStaveRingFixPoint()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLSupport,"STAVERINGFIXINGPOS",index)) {
    double value = db()->getDouble(m_PixelIBLSupport,"STAVERINGFIXINGPOS",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLSupportMidRingWidth()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLSupport,"STAVEMIDRINGWIDTH",index)) {
    double value = db()->getDouble(m_PixelIBLSupport,"STAVEMIDRINGWIDTH",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLSupportMidRingInnerRadius()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLSupport,"STAVEMIDRINGINNERRADIUS",index)) {
    double value = db()->getDouble(m_PixelIBLSupport,"STAVEMIDRINGINNERRADIUS",index) * Gaudi::Units::mm;
    if (value > 0)  return value;
  } 
  return 0.0;
}

double OraclePixGeoManager:: IBLSupportMidRingOuterRadius()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLSupport,"STAVEMIDRINGOUTERRADIUS",index)) {
    double value = db()->getDouble(m_PixelIBLSupport,"STAVEMIDRINGOUTERRADIUS",index) * Gaudi::Units::mm;
    if (value > 0)  return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLFlexMiddleGap()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,"FLEXMIDGAP",index)) {
    double value = db()->getDouble(m_PixelIBLFlex,"FLEXMIDGAP",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

bool OraclePixGeoManager::IBLFlexAndWingDefined()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  return db()->testField(m_PixelIBLFlex,"FLEXMIDGAP",index);
}


double OraclePixGeoManager::IBLFlexDoglegLength()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,"FLEXDOGLEGLENGTH",index)) {
    double value = db()->getDouble(m_PixelIBLFlex,"FLEXDOGLEGLENGTH",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}


double OraclePixGeoManager::IBLStaveFlexWingWidth()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,"FLEXWINGWIDTH",index)) {
    double value = db()->getDouble(m_PixelIBLFlex,"FLEXWINGWIDTH",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLStaveFlexWingThick()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,"FLEXWINGTHICK",index)) {
    double value = db()->getDouble(m_PixelIBLFlex,"FLEXWINGTHICK",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLFlexDoglegRatio()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,"FLEXDOGLEGRATIO",index)) {
    double value = db()->getDouble(m_PixelIBLFlex,"FLEXDOGLEGRATIO",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLFlexDoglegHeight(int iHeight)
{
  std::ostringstream lname;
  lname << "FLEXDOGLEGHEIGHT"<<iHeight;

  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,lname.str(),index)) {
    double value = db()->getDouble(m_PixelIBLFlex,lname.str(),index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLFlexDoglegDY()
{
  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,"FLEXDOGLEGDY",index)) {
    double value = db()->getDouble(m_PixelIBLFlex,"FLEXDOGLEGDY",index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLFlexPP0Z(int iPos)
{
  std::ostringstream lname;
  lname << "FLEXPP0_Z"<<iPos;

  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,lname.str(),index)) {
    double value = db()->getDouble(m_PixelIBLFlex,lname.str(),index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}


double OraclePixGeoManager::IBLFlexPP0Rmin(int iPos)
{
  std::ostringstream lname;
  lname << "FLEXPP0_S"<<iPos<<"RMIN";

  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,lname.str(),index)) {
    double value = db()->getDouble(m_PixelIBLFlex,lname.str(),index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}

double OraclePixGeoManager::IBLFlexPP0Rmax(int iPos)
{
  std::ostringstream lname;
  lname << "FLEXPP0_S"<<iPos<<"RMAX";

  //  int index = PixelStaveIndex(m_currentLD);
  int index=0;
  if (db()->testField(m_PixelIBLFlex,lname.str(),index)) {
    double value = db()->getDouble(m_PixelIBLFlex,lname.str(),index) * Gaudi::Units::mm;
    return value ;
  } 
  return 0.0;
}


std::string OraclePixGeoManager::IBLFlexMaterial(int iPos, const std::string& flexType)
{

  int nbMaterial=db()->getTableSize(m_PixelIBLFlexMaterial);
  int cmptType=0;

  for(int index=0; index<nbMaterial; index++)
    {
      std::string flexTypeIdx = db()->getString(m_PixelIBLFlexMaterial,"TYPE",index);
      if(flexTypeIdx.compare(flexType)==0)
	{
	  cmptType++;
	  if(iPos==cmptType){
	    std::string matTypeIdx = db()->getString(m_PixelIBLFlexMaterial,"MATERIALNAME",index);
	    return matTypeIdx;
	  }
	}
    }
  return std::string("noMat");
      
}


double OraclePixGeoManager:: IBLServiceGetMinRadialPosition(const std::string& srvName, const std::string& srvType, 
							    double srvZmin, double srvZmax)
{
  
  double rmin=99999.;
  
  int numServices =  PixelServiceElements(srvType);
  int nbSrv=0;
  for(int ii = 0; ii < numServices; ii++) {
    // Retrieve/calculate the parameters for the volume.
    //
    std::string name;
    if(srvType=="simple")
      name=db()->getString(m_PixelSimpleService,"NAME",ii);
    else
      name=PixelServiceName(srvType,ii);

    if(name.find(srvName)!=std::string::npos){
      double zmin, zmax, r;
      int symm;
      if(srvType=="simple"){
	zmin=db()->getDouble(m_PixelSimpleService,"ZMIN",ii)*Gaudi::Units::mm;
	zmax=db()->getDouble(m_PixelSimpleService,"ZMAX",ii)*Gaudi::Units::mm;
	symm=db()->getInt(m_PixelSimpleService,"ZSYMM",ii);
	r=db()->getDouble(m_PixelSimpleService,"RMAX",ii)*Gaudi::Units::mm;
      }
      else {
	zmin=PixelServiceZMin(srvType, ii);
	zmax=PixelServiceZMax(srvType, ii);
	symm=PixelServiceZsymm(srvType, ii);
	r=PixelServiceRMin(srvType, ii);
      }

      bool bZintervalle = false;
      if( (srvZmin-zmin)*(srvZmin-zmax)<0 || (srvZmax-zmin)*(srvZmax-zmax)<0 ) bZintervalle=true; 
      if( symm==1 && ((-srvZmin-zmin)*(-srvZmin-zmax)<0 || (-srvZmax-zmin)*(-srvZmax-zmax)<0) ) bZintervalle=true; 
      
      if(bZintervalle){
	if(r<rmin) rmin=r;
	nbSrv++;
      }
    }
  }
  
 if(nbSrv<1)return -1;
  return rmin;

}

double OraclePixGeoManager:: IBLServiceGetMaxRadialPosition(const std::string& srvName, const std::string& srvType, 
							    double srvZmin, double srvZmax)
{
  
  double rmax=-1.;
  int numServices =  PixelServiceElements(srvType);

  int nbSrv=0;
  for(int ii = 0; ii < numServices; ii++) {
    // Retrieve/calculate the parameters for the volume.
    //
    std::string name;
    if(srvType=="simple")
      name=db()->getString(m_PixelSimpleService,"NAME",ii);
    else
      name=PixelServiceName(srvType,ii);

    if(name.find(srvName)!=std::string::npos){

      double zmin, zmax, r;
      int symm;
      if(srvType=="simple"){
	zmin=db()->getDouble(m_PixelSimpleService,"ZMIN",ii)*Gaudi::Units::mm;
	zmax=db()->getDouble(m_PixelSimpleService,"ZMAX",ii)*Gaudi::Units::mm;
	symm=db()->getInt(m_PixelSimpleService,"ZSYMM",ii);
	r=db()->getDouble(m_PixelSimpleService,"RMAX",ii)*Gaudi::Units::mm;
      }
      else {
	zmin=PixelServiceZMin(srvType, ii);
	zmax=PixelServiceZMax(srvType, ii);
	symm=PixelServiceZsymm(srvType, ii);
	r=PixelServiceRMax(srvType, ii);
      }
      
      bool bZintervalle = false;
      if( (srvZmin-zmin)*(srvZmin-zmax)<0 || (srvZmax-zmin)*(srvZmax-zmax)<0 ) bZintervalle=true; 
      if( symm==1 && ((-srvZmin-zmin)*(-srvZmin-zmax)<0 || (-srvZmax-zmin)*(-srvZmax-zmax)<0) ) bZintervalle=true; 
      
      if(bZintervalle && r>rmax){
	rmax=r;
	nbSrv++;
      }
    }
  }
  
 if(nbSrv<1)return -1;
  return rmax;

}

int OraclePixGeoManager::PixelBiStaveType(int layer, int phi)
{
  if (m_staveTypeTable->size() == 0) return phi % 2; 
  return m_pixelStaveTypes->getBiStaveType(layer, phi) % 2;
}

int OraclePixGeoManager::NPixelSectors() 
{
  return db()->getInt(m_PixelLayer,"NSECTORS",m_currentLD);
}

double OraclePixGeoManager::PhiOfModuleZero()
{
  // For backward compatibilty first module is at 1/2 a module division
  if (!db()->testField(m_PixelLayer,"PHIOFMODULEZERO",m_currentLD)){
    if(NPixelSectors()>0) return 180.0*Gaudi::Units::degree/NPixelSectors();
    return 0.;
  } else { 
    return db()->getDouble(m_PixelLayer,"PHIOFMODULEZERO",m_currentLD) * Gaudi::Units::degree;
  }
}


int OraclePixGeoManager::PixelNModule() 
{
  int staveIndex = PixelStaveIndex(m_currentLD);
  if(ibl() && PixelStaveLayout()>3 && PixelStaveLayout()<7 && m_currentLD==0)
    return IBLStaveModuleNumber();
  else
    return db()->getInt(m_PixelStave,"NMODULE",staveIndex);

}

double OraclePixGeoManager::PixelModuleAngle() 
{
  int staveIndex = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"MODULETILT",staveIndex)*Gaudi::Units::deg;
}

double OraclePixGeoManager::PixelModuleDrDistance() 
{
  int staveIndex = PixelStaveIndex(m_currentLD);
  return db()->getDouble(m_PixelStave,"CENTRMODULESHIFT",staveIndex)*mmcm();
}

double OraclePixGeoManager::PixelModuleZPosition(int etaModule) 
{
  // ZPOSTYPE 0. Means equi-distant modules.
  // ZPOSTYPE != 0. Means tabulated z positions.
  int staveIndex = PixelStaveIndex(m_currentLD);
  int zPosType = 0;
  if (ibl() && db()->testField(m_PixelStave,"ZPOSTYPE",staveIndex)) {
    zPosType = db()->getInt(m_PixelStave,"ZPOSTYPE",staveIndex);
  }
  if (zPosType) {
    // Z positions from table
    return PixelModuleZPositionTabulated(etaModule, zPosType);
  } else {
   // Equi-distant modules
    int moduleIndex =  PixelModuleIndexFromEta(etaModule);  
    return db()->getDouble(m_PixelStave,"MODULEDZ",staveIndex)*mmcm() * (moduleIndex - 0.5*(PixelNModule()-1));
  }
}

double OraclePixGeoManager::PixelModuleZPositionTabulated(int etaModule, int type) 
{ 
  if (!m_zPositionMap) {
    m_zPositionMap = new InDetDD::PairIndexMap;
    for (unsigned int indexTmp = 0; indexTmp < db()->getTableSize(m_PixelStaveZ); ++indexTmp) {
      int eta_module = db()->getInt(m_PixelStaveZ,"ETAMODULE",indexTmp);
      int type_tmp       = db()->getInt(m_PixelStaveZ,"TYPE",indexTmp);
      m_zPositionMap->add(type_tmp,eta_module,indexTmp);
    }
  }
  int index = m_zPositionMap->find(type, etaModule);
  if (index < 0) {
    msg(MSG::ERROR) << "Z position not found for etaModule,type =  " << etaModule << ", " << type << endmsg;
    return 0;
  }
  return db()->getDouble(m_PixelStaveZ,"ZPOS",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelModuleShiftFlag(int etaModule) 
{
  if (centerModule(etaModule)) return 1;
  return 0.;
}

double OraclePixGeoManager::PixelModuleStaggerDistance()
{
  int staveIndex = PixelStaveIndex(m_currentLD);
  if (!ibl() || !db()->testField(m_PixelStave,"STAGGERDIST",staveIndex)) return 0; 
  return db()->getDouble(m_PixelStave,"STAGGERDIST",staveIndex) * Gaudi::Units::mm;
}

int OraclePixGeoManager::PixelModuleStaggerSign(int etaModule)
{
  int staveIndex = PixelStaveIndex(m_currentLD);
  if (!ibl() || !db()->testField(m_PixelStave,"FIRSTSTAGGER",staveIndex)) return 0;  
  // FIRSTSTAGGER refers to whether the first module (lowest etavalue) is staggered up (+1) or down(-1)
  int firstStagger =  db()->getInt(m_PixelStave,"FIRSTSTAGGER",staveIndex);
  int moduleIndex = PixelModuleIndexFromEta(etaModule);
  return firstStagger * (moduleIndex%2 ? -1 : 1);
}

bool OraclePixGeoManager::allowSkipEtaZero()
{
  bool allowSkip = true;
  if (ibl()){
    int staveIndex = PixelStaveIndex(m_currentLD);
    if (db()->testField(m_PixelStave,"NOSKIPZERO",staveIndex)) {
      if (db()->getInt(m_PixelStave,"NOSKIPZERO",staveIndex)) allowSkip = false;
    }
  }
  return allowSkip;
}

bool OraclePixGeoManager::centerModule(int etaModule) 
{
  // There is only a center module if there are an odd number
  // of modules. In that case it will be etaModule = 0.
  return (etaModule == 0 && PixelNModule()%2);
}

int OraclePixGeoManager::PixelModuleEtaFromIndex(int index) 
{
  int nModules = PixelNModule();
  int etaModule = index-nModules/2;
  // If even number of modules skip eta = 0.
  // For IBL this behaviour can be disabled.
  if (allowSkipEtaZero() && (etaModule >= 0) && !(nModules%2)) etaModule++; 
  return etaModule;
}

int OraclePixGeoManager::PixelModuleIndexFromEta(int etaModule) 
{
  int nModules = PixelNModule();  
  int index = etaModule + nModules/2;
  // If even number of modules skip eta = 0.
  // For IBL this behaviour can be disabled.
  if (allowSkipEtaZero() && (etaModule >= 0) && (nModules%2 == 0)) index--; 
  return index;
}


double OraclePixGeoManager::PixelModuleAngleSign(int etaModule) 
{
  if (centerModule(etaModule)) return 0;
  if(etaModule < 0) return 1.;
  return -1.;
}

int OraclePixGeoManager::PixelEndcapNDisk() 
{
  return db()->getInt(m_PixelEndcapGeneral,"NDISK");
}

// Endcap container
double  OraclePixGeoManager::PixelEndcapRMin()
{
  return db()->getDouble(m_PixelEndcapGeneral,"RMIN")*mmcm();
}

double  OraclePixGeoManager::PixelEndcapRMax() 
{
  return db()->getDouble(m_PixelEndcapGeneral,"RMAX")*mmcm();
}

double  OraclePixGeoManager::PixelEndcapZMin() 
{
  return db()->getDouble(m_PixelEndcapGeneral,"ZMIN")*mmcm();
}

double  OraclePixGeoManager::PixelEndcapZMax()
{
  return db()->getDouble(m_PixelEndcapGeneral,"ZMAX")*mmcm();
}

int OraclePixGeoManager::PixelEndcapNSupportFrames()
{
   // Obsolete - retus 0 in recent versions
 return (int) db()->getDouble(m_PixelEndcapGeneral,"NFRAME");
}

// Endcap Inner 
double  OraclePixGeoManager::PixelDiskZPosition() 
{
  return db()->getDouble(m_PixelDisk,"ZDISK",m_currentLD)*mmcm();
}

double OraclePixGeoManager::PixelECSiDz1() 
{
  return db()->getDouble(m_PixelDisk,"DZCOUNTER",m_currentLD)*mmcm();
}

double OraclePixGeoManager::PixelECSiDz2() 
{
  return PixelECSiDz1();
}

int OraclePixGeoManager::PixelECNSectors1()
{
  return db()->getInt(m_PixelDisk,"NMODULE",m_currentLD);
}

int OraclePixGeoManager::PixelECNSectors2() 
{
  return PixelECNSectors1();
}

// Endcap Cables
double OraclePixGeoManager::PixelECCablesRMin()
{
  return db()->getDouble(m_PixelDisk,"RMINCABLE",m_currentLD)*mmcm();
}

double OraclePixGeoManager::PixelECCablesRMax()
{
  return db()->getDouble(m_PixelDisk,"RMAXCABLE",m_currentLD)*mmcm();
}


double OraclePixGeoManager::PixelECCablesDistance()
{
  return db()->getDouble(m_PixelDisk,"ZCABLE",m_currentLD)*mmcm();
}

//
/// TMT
//
int OraclePixGeoManager::PixelTMTNumParts()
{
  if (useLegacy()) return m_legacyManager->PixelTMTNumParts();
  return db()->getTableSize(m_PixelTMT);
}

double OraclePixGeoManager::PixelTMTWidthX1(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTWidthX1(iPart);
  return db()->getDouble(m_PixelTMT,"WIDTHX1",iPart) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelTMTWidthX2(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTWidthX2(iPart);
  return db()->getDouble(m_PixelTMT,"WIDTHX2",iPart) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelTMTWidthY(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTWidthY(iPart);
  return db()->getDouble(m_PixelTMT,"WIDTHY",iPart) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelTMTBaseX1(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTBaseX1(iPart);
  return db()->getDouble(m_PixelTMT,"BASEX1",iPart) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelTMTBaseX2(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTBaseX2(iPart);
  return db()->getDouble(m_PixelTMT,"BASEX2",iPart) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelTMTPosY(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTPosY(iPart);
  return db()->getDouble(m_PixelTMT,"Y",iPart) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelTMTPosZ1(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTPosZ1(iPart);
  return db()->getDouble(m_PixelTMT,"Z1",iPart) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelTMTPosZ2(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTPosZ2(iPart);
  return db()->getDouble(m_PixelTMT,"Z2",iPart) * Gaudi::Units::mm;
}

bool OraclePixGeoManager::PixelTMTPerModule(int iPart)
{
  if (useLegacy()) return m_legacyManager->PixelTMTPerModule(iPart);
  return db()->getInt(m_PixelTMT,"PERMODULE",iPart);
}

//
// Omega parameters
//
double OraclePixGeoManager::PixelOmegaUpperBendX()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaUpperBendX();
  return db()->getDouble(m_PixelOmega,"UPPERBENDX") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaUpperBendY()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaUpperBendY();
  return db()->getDouble(m_PixelOmega,"UPPERBENDY") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaUpperBendRadius()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaUpperBendRadius();
  return db()->getDouble(m_PixelOmega,"UPPERBENDR") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaLowerBendX()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaLowerBendX();
  return db()->getDouble(m_PixelOmega,"LOWERBENDX") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaLowerBendY()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaLowerBendY();
  return db()->getDouble(m_PixelOmega,"LOWERBENDY") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaLowerBendRadius()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaLowerBendRadius();
  return db()->getDouble(m_PixelOmega,"LOWERBENDR") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaWallThickness()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaWallThickness();
  return db()->getDouble(m_PixelOmega,"THICK") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaLength()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaLength();
  return db()->getDouble(m_PixelOmega,"LENGTH") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaStartY()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaStartY();
  return db()->getDouble(m_PixelOmega,"STARTY") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaEndY()
{
  if (useLegacy()) return m_legacyManager->PixelOmegaEndY();
  return db()->getDouble(m_PixelOmega,"ENDY") * Gaudi::Units::mm;
}

//
// Al Tube
//

double OraclePixGeoManager::PixelAlTubeUpperBendX()
{
  if (useLegacy()) return m_legacyManager->PixelAlTubeUpperBendX();
  return db()->getDouble(m_PixelAlTube,"UPPERBENDX") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelAlTubeUpperBendY()
{
  if (useLegacy()) return m_legacyManager->PixelAlTubeUpperBendY();
  return db()->getDouble(m_PixelAlTube,"UPPERBENDY") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelAlTubeUpperBendRadius()
{
  if (useLegacy()) return m_legacyManager->PixelAlTubeUpperBendRadius();
  return db()->getDouble(m_PixelAlTube,"UPPERBENDR") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelAlTubeLowerBendX()
{
  if (useLegacy()) return m_legacyManager->PixelAlTubeLowerBendX();
  return db()->getDouble(m_PixelAlTube,"LOWERBENDX") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelAlTubeLowerBendY()
{
  if (useLegacy()) return m_legacyManager->PixelAlTubeLowerBendY();
  return db()->getDouble(m_PixelAlTube,"LOWERBENDY") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelAlTubeLowerBendRadius()
{
  if (useLegacy()) return m_legacyManager->PixelAlTubeLowerBendRadius();
  return db()->getDouble(m_PixelAlTube,"LOWERBENDR") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelAlTubeWallThickness()
{
  if (useLegacy()) return m_legacyManager->PixelAlTubeWallThickness();
  return db()->getDouble(m_PixelAlTube,"THICK") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelAlTubeLength()
{
  if (useLegacy()) return m_legacyManager->PixelAlTubeLength();
  return db()->getDouble(m_PixelAlTube,"LENGTH") * Gaudi::Units::mm;
}

//
// Glue
// 

int OraclePixGeoManager::PixelNumOmegaGlueElements()
{
  if (useLegacy()) return m_legacyManager->PixelNumOmegaGlueElements();
  return db()->getTableSize(m_PixelOmegaGlue);
}

double OraclePixGeoManager::PixelOmegaGlueStartX(int index)
{
  if (useLegacy()) return m_legacyManager->PixelOmegaGlueStartX(index);
  return db()->getDouble(m_PixelOmegaGlue,"STARTX",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaGlueThickness(int index)
{
  if (useLegacy()) return m_legacyManager->PixelOmegaGlueThickness(index);
  return db()->getDouble(m_PixelOmegaGlue,"THICK",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaGlueStartY(int index)
{
  if (useLegacy()) return m_legacyManager->PixelOmegaGlueStartY(index);
  return db()->getDouble(m_PixelOmegaGlue,"STARTY",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaGlueEndY(int index)
{
  if (useLegacy()) return m_legacyManager->PixelOmegaGlueEndY(index);
  return db()->getDouble(m_PixelOmegaGlue,"ENDY",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaGlueLength(int index)
{
  if (useLegacy()) return m_legacyManager->PixelOmegaGlueLength(index);
  return db()->getDouble(m_PixelOmegaGlue,"LENGTH",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelOmegaGluePosZ(int index)
{
  if (useLegacy()) return m_legacyManager->PixelOmegaGluePosZ(index);
  return db()->getDouble(m_PixelOmegaGlue,"Z",index) * Gaudi::Units::mm;
}

int OraclePixGeoManager::PixelOmegaGlueTypeNum(int index)
{
  if (useLegacy()) return m_legacyManager->PixelOmegaGlueTypeNum(index);
  return db()->getInt(m_PixelOmegaGlue,"TYPENUM",index);
}


//
// Fluid
// 
double OraclePixGeoManager::PixelFluidZ1(int index)
{
  if (useLegacy()) return m_legacyManager->PixelFluidZ1(index);
  return db()->getDouble(m_PixelFluid,"Z1",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFluidZ2(int index)
{
  if (useLegacy()) return m_legacyManager->PixelFluidZ2(index);
  return db()->getDouble(m_PixelFluid,"Z2",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFluidThick1(int index)
{
  if (useLegacy()) return m_legacyManager->PixelFluidThick1(index);
  return db()->getDouble(m_PixelFluid,"THICK1",index) * Gaudi::Units::mm;
}


double OraclePixGeoManager::PixelFluidThick2(int index)
{
  if (useLegacy()) return m_legacyManager->PixelFluidThick2(index);
  return db()->getDouble(m_PixelFluid,"THICK2",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFluidWidth(int index)
{
  if (useLegacy()) return m_legacyManager->PixelFluidWidth(index);
  return db()->getDouble(m_PixelFluid,"WIDTH",index) * Gaudi::Units::mm;
}


double OraclePixGeoManager::PixelFluidX(int index)
{
  if (useLegacy()) return m_legacyManager->PixelFluidX(index);
  return db()->getDouble(m_PixelFluid,"X",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelFluidY(int index)
{
  if (useLegacy()) return m_legacyManager->PixelFluidY(index);
  return db()->getDouble(m_PixelFluid,"Y",index) * Gaudi::Units::mm;
}

int OraclePixGeoManager::PixelFluidType(int index)
{
  if (useLegacy()) return m_legacyManager->PixelFluidType(index);
  return db()->getInt(m_PixelFluid,"TYPE",index);
}

int OraclePixGeoManager::PixelFluidNumTypes()
{
  if (useLegacy()) return m_legacyManager->PixelFluidNumTypes();
  return db()->getTableSize(m_PixelFluid);
}

int OraclePixGeoManager::PixelFluidIndex(int type)
{
  for (int i = 0; i < PixelFluidNumTypes(); i++) {
    if (type == PixelFluidType(i)) return i;
  }
  msg(MSG::ERROR) << "Unrecognized fluid volume type: " << type << endmsg;
  return -1;
}
 
std::string OraclePixGeoManager::PixelFluidMat(int index) {
  int matType = 0;
  if (useLegacy()) {
    matType = m_legacyManager->PixelFluidMatType(index);
  } else {
    matType = db()->getInt(m_PixelFluid,"MATTYPE",index);
  }
  return getMaterialName("Fluid", 0, matType);
}

int OraclePixGeoManager::PixelFluidOrient(int layer, int phi) 
{
  if (useLegacy()) return m_legacyManager->PixelFluidOrient(layer, phi);
  return m_pixelStaveTypes->getFluidType(layer,phi);
}

//
// Pigtail
//
double OraclePixGeoManager::PixelPigtailThickness()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailThickness();
  return db()->getDouble(m_PixelPigtail,"THICK") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelPigtailStartY()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailStartY();
  return db()->getDouble(m_PixelPigtail,"STARTY") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelPigtailEndY()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailEndY();
  return db()->getDouble(m_PixelPigtail,"ENDY") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelPigtailWidthZ()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailWidthZ();
  return db()->getDouble(m_PixelPigtail,"WIDTHZ") * Gaudi::Units::mm;
}

// Different width from the curved section in old geometry
double OraclePixGeoManager::PixelPigtailFlatWidthZ()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailFlatWidthZ();
  return PixelPigtailWidthZ();
}

double OraclePixGeoManager::PixelPigtailPosX()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailPosX();
  return db()->getDouble(m_PixelPigtail,"X") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelPigtailPosZ()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailPosZ();
  return db()->getDouble(m_PixelPigtail,"Z") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelPigtailBendX()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailBendX();
  return db()->getDouble(m_PixelPigtail,"BENDX") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelPigtailBendY()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailBendY();
  return db()->getDouble(m_PixelPigtail,"BENDY") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelPigtailBendRMin()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailBendRMin();
  return db()->getDouble(m_PixelPigtail,"BENDRMIN") * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelPigtailBendRMax()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailBendRMax();
  return PixelPigtailBendRMin() + PixelPigtailThickness();
}

double OraclePixGeoManager::PixelPigtailBendPhiMin()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailBendPhiMin();
  return db()->getDouble(m_PixelPigtail,"BENDPHIMIN") * Gaudi::Units::deg;
}

double OraclePixGeoManager::PixelPigtailBendPhiMax()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailBendPhiMax();
  return db()->getDouble(m_PixelPigtail,"BENDPHIMAX") * Gaudi::Units::deg;
}

double OraclePixGeoManager::PixelPigtailEnvelopeLength()
{
  if (useLegacy()) return m_legacyManager->PixelPigtailEnvelopeLength();
  return db()->getDouble(m_PixelPigtail,"ENVLENGTH") * Gaudi::Units::mm;
}

//
// Connector
//
int OraclePixGeoManager::PixelNumConnectorElements()
{
  if (useLegacy()) return m_legacyManager->PixelNumConnectorElements();
  return db()->getTableSize(m_PixelConnector);
}

double OraclePixGeoManager::PixelConnectorWidthX(int index)
{
  if (useLegacy()) return m_legacyManager->PixelConnectorWidthX(index);
  return db()->getDouble(m_PixelConnector,"WIDTHX",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelConnectorWidthY(int index)
{
  if (useLegacy()) return m_legacyManager->PixelConnectorWidthY(index);
  return db()->getDouble(m_PixelConnector,"WIDTHY",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelConnectorWidthZ(int index)
{
  if (useLegacy()) return m_legacyManager->PixelConnectorWidthZ(index);
  return db()->getDouble(m_PixelConnector,"WIDTHZ",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelConnectorPosX(int index)
{
  if (useLegacy()) return m_legacyManager->PixelConnectorPosX(index);
  return db()->getDouble(m_PixelConnector,"X",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelConnectorPosY(int index)
{
  if (useLegacy()) return m_legacyManager->PixelConnectorPosY(index);
  return db()->getDouble(m_PixelConnector,"Y",index) * Gaudi::Units::mm;
}

double OraclePixGeoManager::PixelConnectorPosZ(int index)
{
  if (useLegacy()) return m_legacyManager->PixelConnectorPosZ(index);
  return db()->getDouble(m_PixelConnector,"Z",index) * Gaudi::Units::mm;
}

//
// Design
//


int  OraclePixGeoManager::designType(bool isModule3D)
{
  
  if(isModule3D) return designType3D();

  if (dbVersion() < 2) {
    return 0; // Not used
  } else {
    if (m_PixelReadout->size() == 0) {
      msg(MSG::ERROR) << "ERROR in PixelReadout size. Should not occur!" << endmsg;
      return 0;
    } else if (m_PixelReadout->size() == 1 && !ibl()) {
      return 0;
    } else { // Only in IBL
      return db()->getInt(m_PixelModule,"DESIGNTYPE",moduleType());
    }
  }
}


int  OraclePixGeoManager::designType3D()
{

  if (dbVersion() < 2) {
    return 0; // Not used
  } 
  else {
    if (m_PixelReadout->size() == 0) {
      msg(MSG::ERROR) << "ERROR in PixelReadout size. Should not occur!" << endmsg;
      return 0;
    } else if (m_PixelReadout->size() == 1 && !ibl()) {
      return 0;
    } else { // Only in IBL
      int type = db()->getInt(m_PixelModule,"DESIGNTYPE",moduleType3D());
      return type;
    }
  }
}

int OraclePixGeoManager::DesignReadoutSide(bool isModule3D)
{
  if (dbVersion() < 2) {
    return -1; 
  } else {
    int type = designType((ibl()&&isModule3D));

    return db()->getInt(m_PixelReadout,"READOUTSIDE",type);
  }
}

int OraclePixGeoManager::DesignNumChipsPhi(bool isModule3D)
{
  if (dbVersion() < 2) {  
    return m_legacyManager->DesignNumChipsPhi();
  } else {
    int type = designType((ibl()&&isModule3D));

    return db()->getInt(m_PixelReadout,"NCHIPSPHI",type);
  } 
}    


int OraclePixGeoManager::DesignNumChipsEta(bool isModule3D)
 {
  if (dbVersion() < 2) {  
    return m_legacyManager->DesignNumChipsEta();
  } else {
    int type = designType((ibl()&&isModule3D));

    return db()->getInt(m_PixelReadout,"NCHIPSETA",type);
  }
}

int OraclePixGeoManager::DesignNumRowsPerChip(bool isModule3D)
{
  if (dbVersion() < 2) {  
    return m_legacyManager->DesignNumRowsPerChip(isInnermostPixelLayer());
  } else {
    int type = designType((ibl()&&isModule3D));

    return db()->getInt(m_PixelReadout,"ROWSPERCHIP",type);
  }
}

int OraclePixGeoManager::DesignNumColsPerChip(bool isModule3D) 
{
  if (dbVersion() < 2) {  
    return m_legacyManager->DesignNumColsPerChip(isInnermostPixelLayer());
  } else {
    int type = designType((ibl()&&isModule3D));

    return db()->getInt(m_PixelReadout,"COLSPERCHIP",type);
  }
}


int OraclePixGeoManager::DesignDiodesPhiTotal(bool isModule3D)
{
  if  (dbVersion() < 2) {
    return m_legacyManager->DesignDiodesPhiTotal(isInnermostPixelLayer());
  } else {
    return DesignNumChipsPhi(isModule3D) * (DesignNumRowsPerChip(isModule3D)+DesignNumEmptyRowsInGap(isModule3D)) - DesignNumEmptyRowsInGap(isModule3D);
  }
}

int OraclePixGeoManager::DesignDiodesEtaTotal(bool isModule3D)
{
  if  (dbVersion() < 2) {
    return m_legacyManager->DesignDiodesEtaTotal(isInnermostPixelLayer());
  } else {
    return DesignNumChipsEta(isModule3D) * DesignNumColsPerChip(isModule3D);
  }
}


int OraclePixGeoManager::DesignCellRowsPerCircuit(bool isModule3D)
{
  return DesignNumChipsPhi(isModule3D) * DesignNumRowsPerChip(isModule3D);
}

int OraclePixGeoManager::DesignCellColumnsPerCircuit(bool isModule3D)
{
  return DesignNumColsPerChip(isModule3D);
}

int OraclePixGeoManager::DesignDiodeRowsPerCircuit(bool isModule3D)
{
  return DesignDiodesPhiTotal(isModule3D);
}

int OraclePixGeoManager::DesignDiodeColumnsPerCircuit(bool isModule3D)
{
  return DesignNumColsPerChip(isModule3D);
}

int  OraclePixGeoManager::DesignNumEmptyRowsInGap(bool isModule3D)
{
  // Could determine it from m_gangedIndexMap but expect it to be filled correctly in PixelReadoutTable 
  if (dbVersion() < 2) {
    return m_legacyManager->DesignNumEmptyRowsInGap();
  } else {
    int type=designType((ibl()&&isModule3D));

    return db()->getInt(m_PixelReadout,"EMPTYROWS",type);
  } 
}

// Ganged Pixels
int OraclePixGeoManager::GangedType()
{
  // type 0 means no ganged pixels
  if (!ibl()) return 1;
  if (ibl()) {
    return db()->getInt(m_PixelReadout,"GANGEDTYPE",designType());
  } else {
    int type = 1;
    if (db()->testField(m_PixelReadout,"GANGEDTYPE",designType())) {
      type = db()->getInt(m_PixelReadout,"GANGEDTYPE",designType());
    }
    return type;
  }
}
      

int OraclePixGeoManager::GangedTableIndex(int index, int type)
{
  // There is only one type for standard ATLAS so we just return the index.
  if (!ibl()) return index; 

  if (!m_gangedIndexMap) {
    // First time we create the map
    m_gangedIndexMap = new std::map<int,std::vector<int> >;
    for (unsigned int i = 0; i < db()->getTableSize(m_PixelGangedPixels); i++){
      int testType = 1;
      if (db()->testField(m_PixelGangedPixels,"TYPE",i)) {
	testType = db()->getInt(m_PixelGangedPixels,"TYPE",i);
      }
      (*m_gangedIndexMap)[testType].push_back(i);
    }
  }

  int newIndex = -1;
  std::map<int,std::vector<int> >::const_iterator iter = m_gangedIndexMap->find(type);
  if (iter == m_gangedIndexMap->end()) {
    msg(MSG::ERROR) << "Ganged pixel type " << type << " not found." << endmsg;
  } else {
    const std::vector<int> & vec = iter->second;
    if (index < 0 || static_cast<unsigned int>(index) >= vec.size()) {
      msg(MSG::ERROR) << "Ganged pixel index " << index << " for type " << type << " out of range." << endmsg;
    } else {
      newIndex = vec[index];
    }
  }
  return newIndex;
}
 
int OraclePixGeoManager::NumberOfEmptyRows(bool isModule3D)
{
  return DesignNumEmptyRowsInGap(isModule3D);
}

int OraclePixGeoManager::EmptyRows(int index)
{
  if (dbVersion() < 2) {
    return m_legacyManager->EmptyRows(index);
  } else {
    int newIndex = GangedTableIndex(index, GangedType());
    if (newIndex >= 0) {
      return db()->getInt(m_PixelGangedPixels,"EMPTYROW",newIndex);
    } else {
      return 0;
    }
  }
}

int OraclePixGeoManager::EmptyRowConnections(int index)
{
  if (dbVersion() < 2) {
    return m_legacyManager->EmptyRowConnections(index);
  } else {
    int newIndex = GangedTableIndex(index, GangedType());
    if (newIndex >= 0) {
      return db()->getInt(m_PixelGangedPixels,"CONNECTROW",newIndex);
    } else {
      return 0;
    }
  }    
}


double OraclePixGeoManager::DesignRPActiveArea(bool isModule3D)
{
  if (dbVersion() < 2) { 
    return m_legacyManager->DesignRPActiveArea();
  } else {
    // All layers assumed to be the same.
    return DesignPitchRP(isModule3D) * ((DesignNumRowsPerChip(isModule3D)+DesignNumEmptyRowsInGap(isModule3D)) * DesignNumChipsPhi(isModule3D) - DesignNumEmptyRowsInGap(isModule3D));
  } 
}

double OraclePixGeoManager::DesignZActiveArea(bool isModule3D)
{
  if (dbVersion() < 2) {   
    return m_legacyManager->DesignZActiveArea();
  } else {
    // All layers assumed to be the same.
    return (DesignPitchZ(isModule3D) * (DesignNumColsPerChip(isModule3D) - 2) + 2 * DesignPitchZLong(isModule3D)) * DesignNumChipsEta(isModule3D) + 2 * (DesignPitchZLongEnd(isModule3D) - DesignPitchZLong(isModule3D));
  }
}
    
double OraclePixGeoManager::DesignPitchRP(bool isModule3D)
{
  if (dbVersion() < 2) {
    return m_legacyManager->DesignPitchRP(isInnermostPixelLayer());
  } else {
    int type = designType((ibl()&&isModule3D));
    return db()->getDouble(m_PixelReadout,"PITCHPHI",type) * Gaudi::Units::mm;
 } 
}

double OraclePixGeoManager::DesignPitchZ(bool isModule3D)
{
  if (dbVersion() < 2) {
    return m_legacyManager->DesignPitchZ(isInnermostPixelLayer());
  } else {
    int type = designType((ibl()&&isModule3D));
    return db()->getDouble(m_PixelReadout,"PITCHETA",type) * Gaudi::Units::mm;
  }
}

double OraclePixGeoManager::DesignPitchZLong(bool isModule3D)
{
  // Defaults to DesignPitchZ if not specified or is zero.
  if (dbVersion() < 2) {
    return m_legacyManager->DesignPitchZLong(isInnermostPixelLayer());
  } else {
    int type = designType((ibl()&&isModule3D));
    double pitch = db()->getDouble(m_PixelReadout,"PITCHETALONG",type) * Gaudi::Units::mm;
    if (pitch == 0) pitch = DesignPitchZ(isModule3D);
    return pitch;
  }
}

double OraclePixGeoManager::DesignPitchZLongEnd(bool isModule3D)
{
  // Defaults to DesignPitchZLongEnd if not specified or is zero.
  if (!ibl()) { // This check is not really needed once the field is in the database.
    return DesignPitchZLong(isModule3D);
  } else {
    int type = designType((ibl()&&isModule3D));
    double pitch = 0;
    if (db()->testField(m_PixelReadout,"PITCHETAEND",type)) {
      pitch = db()->getDouble(m_PixelReadout,"PITCHETAEND",type) * Gaudi::Units::mm;
    }
    if (pitch == 0) pitch = DesignPitchZLong(isModule3D);
    return pitch;
  }
}


double OraclePixGeoManager::DesignGapRP(bool isModule3D)
{
  if (dbVersion() < 2) {
    return m_legacyManager->DesignGapRP();
  } else {
    return DesignNumEmptyRowsInGap(isModule3D) * DesignPitchRP(isModule3D);
  }
}

double OraclePixGeoManager::DesignGapZ(bool isModule3D)
{
  if (dbVersion() < 2) {  
    return m_legacyManager->DesignGapZ();
  } else {
    return 2. * (DesignPitchZLong(isModule3D) - DesignPitchZ(isModule3D));
  }
}

int OraclePixGeoManager::DesignCircuitsPhi(bool /* isModule3D */)
{
  //
  // This should be (*pdch)[0]->getDouble("NRPCHIP"), but in the current
  // design we prefer to have one chip in the rphi direction
  // and define the connections for the pixels in the gap
  return 1;
}

int OraclePixGeoManager::DesignCircuitsEta(bool isModule3D)
{
  return DesignNumChipsEta(isModule3D);
}



// Endcap 
double  OraclePixGeoManager::PixelDiskRMin()
{
  return db()->getDouble(m_PixelDisk,"RIDISK",m_currentLD)*mmcm();
}

///

//
// endcap rings
//
int OraclePixGeoManager::PixelDiskNumSupports() {
  // Hardwire for now
  return 3;
}

double OraclePixGeoManager::PixelDiskSupportRMin(int isup) {
  std::ostringstream field;
  field <<"SUP"<< isup+1 <<"RMIN";
  return db()->getDouble(m_PixelDisk,field.str(),m_currentLD)*mmcm();
}

double OraclePixGeoManager::PixelDiskSupportRMax(int isup) {
  std::ostringstream field;
  field <<"SUP"<< isup+1 <<"RMAX";
  return db()->getDouble(m_PixelDisk,field.str(),m_currentLD)*mmcm();
}


// SLHC only (TODO: does not look like it)
double OraclePixGeoManager::PixelDiskSupportThickness(int isup) {

  std::ostringstream prefix;
  prefix <<"SUP"<< isup+1 <<"THICK";

  bool found = false;
  double tck = 0;

  // First check text file
  // default support thickness
  if (db()->testFieldTxt(m_PixelDisk,"SUP_THICK")) {
    tck = db()->getDouble(m_PixelDisk,"SUP_THICK");
    found = true;
  } 
  // overwrites if found
  if (db()->testFieldTxt(m_PixelDisk,prefix.str(),m_currentLD)) {
    tck = db()->getDouble(m_PixelDisk,prefix.str(),m_currentLD);
    found = true;
  }

  // Now check database
  if (!found)  tck = db()->getDouble(m_PixelDisk,prefix.str(),m_currentLD);

  if(tck>0.) {
    return tck * mmcm();
  } else { // radlen
    int typeNum = PixelDiskSupportMaterialTypeNum(isup);
    std::string matName = getMaterialName("DiskSupport", m_currentLD, typeNum);
    return CalculateThickness(tck, matName);
  }
}

// SLHC only (TODO: does not look like it)
int OraclePixGeoManager::PixelDiskSupportMaterialTypeNum(int isup) {
 
  if (dbVersion() < 3) return 0;

  std::ostringstream prefix;
  prefix <<"SUP"<< isup+1 <<"MAT";

  int imat = 0;
  bool found = false;
  // default material type
  if (db()->testFieldTxt(m_PixelDisk,"SUP_MAT")) {
    imat = db()->getInt(m_PixelDisk,"SUP_MAT");
    found = true;
  } 
  // overwrites if found
  if (db()->testFieldTxt(m_PixelDisk,prefix.str(),m_currentLD)) {
    imat = db()->getInt(m_PixelDisk,prefix.str(),m_currentLD);
    found = true;
  }

  if (!found) {
    imat = db()->getInt(m_PixelDisk,prefix.str(),m_currentLD);
  }
  return imat;
}


//
//*** DBM Parameters with local database  ***//
//

// return angle of the telescope
double OraclePixGeoManager::DBMAngle() {
  return db()->getDouble(m_DBMTelescope,"ANGLE")*Gaudi::Units::deg;
}

// return dimension of the DBM telescope
double OraclePixGeoManager::DBMTelescopeX() {
   return db()->getDouble(m_DBMTelescope,"WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMTelescopeY() {
   return db()->getDouble(m_DBMTelescope,"HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMTelescopeZ() {
   return db()->getDouble(m_DBMTelescope,"LENGTH")*Gaudi::Units::mm;
}

// return height and length of the module cage having a 3-layers structure
double OraclePixGeoManager::DBMModuleCageY() {
  return db()->getDouble(m_DBMTelescope,"CAGE_HEIGHT")*Gaudi::Units::mm;
} 
double OraclePixGeoManager::DBMModuleCageZ() {
  return db()->getDouble(m_DBMTelescope,"CAGE_LENGTH")*Gaudi::Units::mm;
} 

// return layer spacing
double OraclePixGeoManager::DBMSpacingZ() {
  return db()->getDouble(m_DBMCage,"ZSPACING")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMSpacingRadial() {
  if (m_currentLD == 0)
    return db()->getDouble(m_DBMCage,"RADIAL_SPACE_0")*Gaudi::Units::mm;
  else if (m_currentLD == 1)
    return db()->getDouble(m_DBMCage,"RADIAL_SPACE_1")*Gaudi::Units::mm;
  else if (m_currentLD == 2)
    return db()->getDouble(m_DBMCage,"RADIAL_SPACE_2")*Gaudi::Units::mm;
  else {
     msg(MSG::WARNING) << "DBMSpacingRadial() is not found" << endmsg;
     return 0.;
  }
}
// return dimension of bracket unit
double OraclePixGeoManager::DBMBracketX() {
  return db()->getDouble(m_DBMBracket,"WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBracketY() {
  return db()->getDouble(m_DBMBracket,"HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBracketZ() {
  return db()->getDouble(m_DBMBracket,"THICKNESS")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMTrapezBackTheta() {
  return db()->getDouble(m_DBMBracket,"TRAPEZBACK_THETA")*Gaudi::Units::deg;
}
double OraclePixGeoManager::DBMTrapezBackX() {
  return db()->getDouble(m_DBMBracket,"TRAPEZBACK_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMTrapezBackY() {
  return db()->getDouble(m_DBMBracket,"TRAPEZBACK_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMTrapezBackShortZ() {
  return db()->getDouble(m_DBMBracket,"TRAPEZBACK_ZSHORT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktWindowX() {
  return db()->getDouble(m_DBMBracket,"WINDOW_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktWindowY() {
  return db()->getDouble(m_DBMBracket,"WINDOW_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktWindowOffset() {
  return db()->getDouble(m_DBMBracket,"WINDOW_OFFSET")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktWindowCenterZ() {
  return db()->getDouble(m_DBMBracket,"WINDOW_CENTERZ")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktTopBlockZ() {
  return db()->getDouble(m_DBMBracket,"TOPBLOCK_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktSideBlockX() {
  return db()->getDouble(m_DBMBracket,"SIDEBLOCK_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktSideBlockY() {
  return db()->getDouble(m_DBMBracket,"SIDEBLOCK_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktLockZ() {
  return db()->getDouble(m_DBMBracket,"LOCK_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktLockY() {
  return db()->getDouble(m_DBMBracket,"LOCK_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktFinLongZ() {
  return db()->getDouble(m_DBMBracket,"COOLINGFIN_ZLONG")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktFinHeight() {
  return db()->getDouble(m_DBMBracket,"COOLINGFIN_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktFinThick() {
  return db()->getDouble(m_DBMBracket,"COOLINGFIN_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMBrcktFinPos() {
  return db()->getDouble(m_DBMBracket,"COOLINGFIN_POS")*Gaudi::Units::mm;
}

// return spacing between V-slide and first layer
double OraclePixGeoManager::DBMSpace() {
  return db()->getDouble(m_DBMCage,"SPACING1")*Gaudi::Units::mm;
}

// return dimensions of the main plate
double OraclePixGeoManager::DBMMainPlateX() {
  return db()->getDouble(m_DBMCage,"MAINPLATE_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMMainPlateY() {
  return db()->getDouble(m_DBMCage,"MAINPLATE_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMMainPlateZ() {
  return db()->getDouble(m_DBMCage,"MAINPLATE_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMMPlateWindowWidth() {
  return db()->getDouble(m_DBMCage,"MPWINDOW_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMMPlateWindowHeight() {
  return db()->getDouble(m_DBMCage,"MPWINDOW_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMMPlateWindowPos() {
  return db()->getDouble(m_DBMCage,"MPWINDOW_POS")*Gaudi::Units::mm;
}
// return dimensions of aluminium side plates
double OraclePixGeoManager::DBMCoolingSidePlateX() {
  return db()->getDouble(m_DBMCage,"SIDEPLATE_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMCoolingSidePlateY() {
  return db()->getDouble(m_DBMCage,"SIDEPLATE_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMCoolingSidePlateZ() {
  return db()->getDouble(m_DBMCage,"SIDEPLATE_LENGTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMCoolingSidePlatePos() {
  return db()->getDouble(m_DBMCage,"SIDEPLATE_POS")*Gaudi::Units::mm;
}

// return dimension of sensor, chip and ceramic
double OraclePixGeoManager::DBMDiamondX() {
  return db()->getDouble(m_DBMModule,"DIAMOND_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMDiamondY() {
  return db()->getDouble(m_DBMModule,"DIAMOND_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMDiamondZ() {
  return db()->getDouble(m_DBMModule,"DIAMOND_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMFEI4X() {
  return db()->getDouble(m_DBMModule,"FEI4_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMFEI4Y() {
  return db()->getDouble(m_DBMModule,"FEI4_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMFEI4Z() {
  return db()->getDouble(m_DBMModule,"FEI4_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMCeramicX() {
  return db()->getDouble(m_DBMModule,"CERAMIC_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMCeramicY() {
  return db()->getDouble(m_DBMModule,"CERAMIC_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMCeramicZ() {
  return db()->getDouble(m_DBMModule,"CERAMIC_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMAirGap() {
  return db()->getDouble(m_DBMModule,"AIR_GAP")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMKaptonZ() {
  return db()->getDouble(m_DBMModule,"KAPTONZ")*Gaudi::Units::mm;
}

// flex support
double OraclePixGeoManager::DBMFlexSupportX() {
  return db()->getDouble(m_DBMCage,"FLEXSUPP_WIDTH")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMFlexSupportY() {
    return db()->getDouble(m_DBMCage,"FLEXSUPP_HEIGHT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMFlexSupportZ() {
  return db()->getDouble(m_DBMCage,"FLEXSUPP_THICK")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMFlexSupportOffset() {
    return db()->getDouble(m_DBMCage, "FLEXSUPP_OFFSET")*Gaudi::Units::mm;
}

// return radius of supporting rod
double OraclePixGeoManager::DBMRodRadius() {
  return db()->getDouble(m_DBMCage,"ROD_RADIUS")*Gaudi::Units::mm;
}
// return distance between center of rods
double OraclePixGeoManager::DBMMPlateRod2RodY() {
  return db()->getDouble(m_DBMCage,"ROD2ROD_VERT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMMPlateRod2RodX() {
  return db()->getDouble(m_DBMCage,"ROD2ROD_HOR")*Gaudi::Units::mm;
}

// radius and thickness of PP0 board
double OraclePixGeoManager::DBMPP0RIn() {
  return db()->getDouble(m_DBMTelescope,"PP0_RIN")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMPP0ROut() {
  return db()->getDouble(m_DBMTelescope,"PP0_ROUT")*Gaudi::Units::mm;
}
double OraclePixGeoManager::DBMPP0Thick() {
  return db()->getDouble(m_DBMTelescope,"PP0_THICK")*Gaudi::Units::mm;
}


