/*
Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/ 
#include "PixelServicesTool/SimpleServiceVolumeMaker.h"
#include "PixelServicesTool/PixelSimpleServiceXMLHelper.h"

#include "InDetGeoModelUtils/ServiceVolume.h"
/*#include "InDetGeoModelUtils/InDetDDAthenaComps.h"*/

#include "InDetGeoModelUtils/OraclePixGeoAccessor.h"

#include "RDBAccessSvc/IRDBRecordset.h"
#include "GeometryDBSvc/IGeometryDBSvc.h"

#include "PathResolver/PathResolver.h"
#include "CLHEP/Units/SystemOfUnits.h"

namespace InDetDD {

SimpleServiceVolumeSchema::SimpleServiceVolumeSchema() 
{
  m_rmin = "RIN";
  m_rmax = "ROUT";
  m_rmin2 = "RIN2";
  m_rmax2 = "ROUT2";
  m_zmin = "ZIN";
  m_zmax = "ZOUT";
  m_zsymm = "ZSYMM";
  m_materialName = "MATERIALNAME";
  m_repeat = "REPEAT";
  m_phiStart = "PHI";
  m_phiDelta = "WIDTH";
  m_width = "WIDTH";
  m_shapeType = "SHAPE";
  m_volName = "VOLNAME";
  m_radialDiv = "";
  m_phiStep = "";
}

SimpleServiceVolumeMakerMgr::SimpleServiceVolumeMakerMgr(const std::string & table, bool readDataFromDB, const PixelGeoBuilderBasics* basics)
  : GeoXMLUtils(),
    PixelGeoBuilder(basics),
    m_table(table),
    m_simpleSrvXMLHelper(0),
    m_readFromDB(readDataFromDB),
    m_XMLdefined(false)
{
  if(!m_readFromDB){
    m_simpleSrvXMLHelper = new PixelSimpleServiceXMLHelper(table, basics);
    m_XMLdefined = true;
  }
}

// Question - why does the db() getter never return the database?
const IGeometryDBSvc *SimpleServiceVolumeMakerMgr::db() const {
  //  return m_athenaComps->geomDB();
  return 0;
}

double SimpleServiceVolumeMakerMgr::rmin(int index) const
{  
  return m_simpleSrvXMLHelper->rmin(index);
}


double SimpleServiceVolumeMakerMgr::rmax(int index) const
{ 
  return m_simpleSrvXMLHelper->rmax(index);
}


double SimpleServiceVolumeMakerMgr::rmin2(int index) const
{
  return m_simpleSrvXMLHelper->rmin2(index);
}

double SimpleServiceVolumeMakerMgr::rmax2(int index) const
{
  return m_simpleSrvXMLHelper->rmax2(index);
}

double SimpleServiceVolumeMakerMgr::zmin(int index) const
{ 
  return m_simpleSrvXMLHelper->zmin(index);
}

double SimpleServiceVolumeMakerMgr::zmax(int index) const
{
  
  return m_simpleSrvXMLHelper->zmax(index);
}

double SimpleServiceVolumeMakerMgr::phiDelta(int index) const
{
  return m_simpleSrvXMLHelper->phiDelta(index);
}

double SimpleServiceVolumeMakerMgr::width(int index) const
{
  return m_simpleSrvXMLHelper->width(index);
}

double SimpleServiceVolumeMakerMgr::phiStart(int index) const
{
  return m_simpleSrvXMLHelper->phiStart(index);
}

double SimpleServiceVolumeMakerMgr::phiStep(int index) const
{
  return 0;
}

bool SimpleServiceVolumeMakerMgr::zsymm(int index) const
{
  return m_simpleSrvXMLHelper->zsymm(index);
}


int SimpleServiceVolumeMakerMgr::repeat(int index) const
{
  return m_simpleSrvXMLHelper->repeat(index);
}

int SimpleServiceVolumeMakerMgr::radialDiv(int index) const
{
  return 0;
}

std::string SimpleServiceVolumeMakerMgr::shapeType(int index) const
{
  return m_simpleSrvXMLHelper->shapeType(index);
}

std::string SimpleServiceVolumeMakerMgr::volName(int index) const
{
  return m_simpleSrvXMLHelper->volName(index);
}

std::string SimpleServiceVolumeMakerMgr::materialName(int index) const
{
  return m_simpleSrvXMLHelper->materialName(index);
}

unsigned int SimpleServiceVolumeMakerMgr::numElements() const
{
  return m_simpleSrvXMLHelper->numElements();
}


SimpleServiceVolumeMaker::SimpleServiceVolumeMaker(const std::string & table, const std::string & label, const PixelGeoBuilderBasics* basics, bool readDataFromDB) 
  : m_table(table),
    m_label(label)
{
  m_mgr = new SimpleServiceVolumeMakerMgr(table, readDataFromDB, basics);
}

SimpleServiceVolumeMaker::~SimpleServiceVolumeMaker()
{
  for (unsigned int i = 0; i < m_services.size(); i++) {
    delete m_services[i];
  }
  delete m_mgr;
}

const std::vector<const ServiceVolume *>& SimpleServiceVolumeMaker::makeAll()
{
  for (unsigned int ii = 0; ii < numElements(); ++ii) {
    m_services.push_back(make(ii));
  }
  return m_services;
}

unsigned int SimpleServiceVolumeMaker::numElements() const {
  return m_mgr->numElements();
}

InDetDD::ServiceVolume *SimpleServiceVolumeMaker::make(int ii)
{
  //
  // Retrieve/calculate the parameters for the volume.
  //
  ServiceVolume * param = new ServiceVolume ;
  param->setMaterial(m_mgr->materialName(ii));
  param->setRmin(m_mgr->rmin(ii));
  param->setRmax(m_mgr->rmax(ii));
  param->setZmin(m_mgr->zmin(ii));
  param->setZmax(m_mgr->zmax(ii));
  param->setZsymm(m_mgr->zsymm(ii));
  param->setVolName(m_mgr->volName(ii));
  
  
  bool needsRotation = false;
  
  // For TUBE there is no need to read the rest 
  std::string shapeType = m_mgr->shapeType(ii);
  if (!shapeType.empty() && shapeType != "TUBE") {
      
    double rmin2 = m_mgr->rmin2(ii);
    double rmax2 = m_mgr->rmax2(ii);
    
    if (rmin2 <= 0) rmin2 = param->rmin(); 
    if (rmax2 <= 0) rmax2 = param->rmax(); 
    
    int radialDiv = m_mgr->radialDiv(ii);
    
    double phiDelta =  m_mgr->phiDelta(ii);
    
    bool fullPhiSector = false;
    if (phiDelta == 0 || phiDelta >=359.9*CLHEP::degree) {
      phiDelta = 360*CLHEP::degree;
      fullPhiSector = true;
    } 
    //else {
    //phiDelta -= 2*phiepsilon;
    //phiStart += phiepsilon;
    // }
    
    if (shapeType == "UNKNOWN") {
      if (radialDiv > 0) {
	shapeType = "RADIAL";
      } else if (param->rmin() == rmin2  &&  param->rmax() == rmax2 ) {
	if (fullPhiSector) {
	  shapeType = "TUBE";
	} else {
	  shapeType = "TUBS";
	}
      } else {
	shapeType = "CONS";
      } 
    }
    
    
    int repeat = m_mgr->repeat(ii);
    if (repeat == 0) repeat = 1;
    
    double phiStart =  m_mgr->phiStart(ii);
    double phiWidth =  phiDelta;
    
    if (shapeType == "CONS"  || shapeType == "TUBS") { 
      const double phiepsilon = 0.001*CLHEP::degree;
      phiWidth -= 2*phiepsilon;
      phiStart += phiepsilon;
    }
    
    // Can be in degree or CLHEP::mm. Usually it is CLHEP::deg expect for BOX, TRAP and ROD shape
    // Geometry manager makes no assumptions about units. So we must interpret here.
    if (shapeType == "BOX" || shapeType == "ROD" || shapeType=="ROD2" || shapeType == "TRAP") {
      phiWidth = m_mgr->width(ii); // in mm
    } 
    
    if (shapeType == "PGON"  || shapeType == "PGON2" || 
	shapeType == "CONE"  || shapeType == "CONS" || 
	shapeType == "PGON3" || shapeType == "PGON4") {
      if ((rmin2 != param->rmin()) || (rmax2 != param->rmax())) {
	needsRotation = true;
      }
    }
    
    int sides = 0;
    int nCopies = 1;
    if (shapeType == "PGON"  || shapeType == "PGON2" ||
	shapeType == "PGON3" || shapeType == "PGON4") {
      sides = repeat;
    } else {
      nCopies = repeat;
    }
    
    // Force nCopies to 1 for TUBE and CONE 
    if (shapeType.empty() || shapeType == "TUBE" || shapeType == "CONE") {
      nCopies = 1;
    }
    
    param->setShapeType(shapeType);
    param->setRmin2(rmin2);
    param->setRmax2(rmax2);
    param->setPhiLoc(phiStart);
    param->setPhiWidth(phiWidth);
    param->setSides(sides);
    param->setNCopies(nCopies);
    //param->setRadialDiv(radialDiv);
    //param->setPhiStep(phiStep);
  }
  
  param->setNeedsRotation(needsRotation);
  
  
  //
  // If zin is 0... (within 10^-5) this is a volume symmetric around
  // the origin
  //
  if(std::abs(param->zmin()) < 0.000001) {
    param->setZmin(-param->zmax());
    param->setZsymm(false);
  }	
  
  param->setLabel(m_label); 

//   std::cout<<"SimpleServiceVoluemMaker "<<ii<<std::endl;
//   param->print();
  
  return param;
}

} // end namespace
