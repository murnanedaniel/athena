/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "RDBMaterialManager.h"
#include "GeoModelUtilities/DecodeVersionKey.h"
#include "GeoModelInterfaces/IGeoModelSvc.h"

#include "GeoModelKernel/GeoMaterial.h"
#include "GeoModelKernel/GeoElement.h"
#include "GeoModelKernel/Units.h"

#include "StoreGate/DataHandle.h"

#include "RDBAccessSvc/IRDBAccessSvc.h"
#include "RDBAccessSvc/IRDBRecordset.h"
#include "RDBAccessSvc/IRDBRecord.h"

#include "AthenaKernel/getMessageSvc.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "AthenaBaseComps/AthCheckMacros.h"
#include "boost/algorithm/string/predicate.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>


//---------------------------Help find elements in the list-----------------------//
class NameEquals {                                                                //
public:                                                                           //
  NameEquals(const std::string & name):m_name(name){}                              //
  bool operator() (const GeoElement *e) const {return m_name==e->getName();}       //
private:                                                                          //
  std::string m_name;                                                              //
};                                                                                //
//--------------------------------------------------------------------------------//

//---------------------------Help find elements in the list-----------------------//
class NumberEquals {                                                              //
public:                                                                           //
  NumberEquals(unsigned int number):m_number(number){}                             //
  bool operator() (const GeoElement *e) const {return m_number==e->getZ();}        //
private:                                                                          //
  unsigned int m_number;                                                           //
};                                                                                //
//--------------------------------------------------------------------------------//

int CheckElement(std::string &name)
{
  if(name.find("::",0) == std::string::npos) {
    return 1;
  }
  else {
    return 0;	
  }
}

int printElement ( GeoElement* &p_element)
{
  std::string name = p_element->getName();
  std::string symbol = p_element->getSymbol();
  double a = p_element->getA();
  double z = p_element->getZ();
	
  std::cout << " ***** CheckElement(): Print the Element:  " << name << std::endl; 
  std::cout << " ***** The Element: name,		symbol, 	A, 	Z " << std::endl; 
  std::cout << " *****             "<<name <<"		"<<symbol <<"		"<< a * (Gaudi::Units::mole / GeoModelKernelUnits::gram) <<"	"<< z <<"	"  << std::endl;
	
  return 1;
}

int printElement ( const GeoElement* &p_element)
{
  std::string name = p_element->getName();
  std::string symbol = p_element->getSymbol();
  double a = p_element->getA();
  double z = p_element->getZ();
	
  std::cout << " ***** PrintElement(): Print the Element:  " << name << std::endl; 
  std::cout << " ***** The Element: name,		symbol, 	A, 	Z " << std::endl; 
  std::cout << " *****             "<<name <<"		"<<symbol <<"		"<< a * (Gaudi::Units::mole / GeoModelKernelUnits::gram) <<"	"<< z <<"	"  << std::endl;
	
  return 1;
}

int printMaterial ( GeoMaterial* &p_material)
{
  std::string name = p_material->getName();
  double density = p_material->getDensity() * (Gaudi::Units::cm3 / GeoModelKernelUnits::gram);

  std::cout << " ***** PrintMaterial(): Print the Material:  " << name << std::endl; 
  std::cout << " ***** The Material: name,	density	" << std::endl; 
  std::cout << " *****              "<< name <<"		"<<density <<"		" << std::endl; 	
	
  return 1;
}

int printFullMaterial ( GeoMaterial* &p_material)
{
  std::string name = p_material->getName();
  double density = p_material->getDensity() * (Gaudi::Units::cm3 / GeoModelKernelUnits::gram);
	
  std::cout << " ***** PrintFullMaterial(): Print the Material:  " << name << std::endl; 
  std::cout << " ***** The Material: name, 	density" << std::endl; 
  std::cout << " *****              "<< name <<" 	 "<<density <<"  " << std::endl; 
	
  p_material->lock();
  int element_number = p_material->getNumElements();	
		
 			
  if ( element_number  == 0){
    std::cout << " ***** No Elements now in this printMaterial( ) " << std::endl;	
    return 1;
  }
  else {
    element_number = p_material->getNumElements();	
	
    for(int i =0; i< element_number;i ++)
      {
	const GeoElement* tmp_element = p_material->getElement(i);
	double element_fraction = p_material->getFraction(i);
		
	std::cout<<" ***** ***** Number:  " << i << " Fraction:  " << element_fraction<< std::endl;
	printElement( tmp_element); 
      }	
  }
  return 1;
}
	
	

RDBMaterialManager::RDBMaterialManager(ISvcLocator* pSvcLocator)
{
  if(!readMaterialsFromDB(pSvcLocator).isSuccess()) {
    throw std::runtime_error("RDBMaterialManager failed to read Geometry DB");
  }
}

StatusCode RDBMaterialManager::readMaterialsFromDB(ISvcLocator* pSvcLocator)
{
  IGeoModelSvc*  iGeoModel;		
  IRDBAccessSvc* iAccessSvc;
  MsgStream log(Athena::getMessageSvc(), "GeoModelSvc::RDBMaterialManager"); 		

  ATH_CHECK(pSvcLocator->service("GeoModelSvc",iGeoModel));
  ATH_CHECK(pSvcLocator->service("RDBAccessSvc",iAccessSvc));

  // Do not load defaults for RUN4
  bool loadDefaults = iGeoModel->geoConfig() != GeoModel::GEO_RUN4;
  log << MSG::DEBUG << "Will load material defaults if not present: " << loadDefaults << endmsg;

  // --- Standard materials, elements
  DecodeVersionKey keyAtlas(iGeoModel, "ATLAS");
  m_elements = iAccessSvc->getRecordsetPtr("Elements",keyAtlas.tag(),keyAtlas.node());
  if(loadDefaults && m_elements->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting Elements with default tag" <<endmsg;
    m_elements = iAccessSvc->getRecordsetPtr("Elements","Materials-00","Materials");
  }
  m_stdmatcomponents = iAccessSvc->getRecordsetPtr("StdMatComponents",keyAtlas.tag(),keyAtlas.node());
  if(loadDefaults && m_stdmatcomponents->size()==0)	{
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting StdMatComponents with default tag" <<endmsg;
    m_stdmatcomponents = iAccessSvc->getRecordsetPtr("StdMatComponents","Materials-00","Materials");
  }
  m_stdmaterials = iAccessSvc->getRecordsetPtr("StdMaterials",keyAtlas.tag(),keyAtlas.node());
  if(loadDefaults && m_stdmaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting StdMaterials with default tag" <<endmsg;
    m_stdmaterials = iAccessSvc->getRecordsetPtr("StdMaterials","Materials-00","Materials");
  }
  
  // --- Pixel materials
  DecodeVersionKey keyPixel(iGeoModel, "Pixel");
  m_pixmatcomponents = iAccessSvc->getRecordsetPtr("PixMatComponents",keyPixel.tag(),keyPixel.node());
  if(loadDefaults && m_pixmatcomponents->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting PixMatComponents with default tag" <<endmsg;
    m_pixmatcomponents = iAccessSvc->getRecordsetPtr("PixMatComponents","PixMatComponents-00");
  }
  m_pixmaterials = iAccessSvc->getRecordsetPtr("PixMaterials",keyPixel.tag(),keyPixel.node());
  if(loadDefaults && m_pixmaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting PixMaterials with default tag" <<endmsg;
    m_pixmaterials = iAccessSvc->getRecordsetPtr("PixMaterials","PixMaterials-00");
  }
  
  // --- Pixel materials for TB
  m_pixtbmatcomponents = iAccessSvc->getRecordsetPtr("PixelTBMatComponents",keyPixel.tag(),keyPixel.node());
  if(loadDefaults && m_pixtbmatcomponents->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting PixTBMatComponents with default tag" <<endmsg;
    m_pixtbmatcomponents = iAccessSvc->getRecordsetPtr("PixMatComponents","PixMatComponents-00");
  }
  m_pixtbmaterials = iAccessSvc->getRecordsetPtr("PixelTBMaterials",keyPixel.tag(),keyPixel.node());
  if(loadDefaults && m_pixtbmaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting PixTBMaterials with default tag" <<endmsg;
    m_pixtbmaterials = iAccessSvc->getRecordsetPtr("PixMaterials","PixMaterials-00");
  }
  
  // --- SCT materials
  DecodeVersionKey keySCT(iGeoModel, "SCT");
  m_sctmatcomponents = iAccessSvc->getRecordsetPtr("SCTMatComponents",keySCT.tag(),keySCT.node());
  if(loadDefaults && m_sctmatcomponents->size()==0)	{
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting SCTMatComponents with default tag" <<endmsg;
    m_sctmatcomponents = iAccessSvc->getRecordsetPtr("SCTMatComponents","SCTMatComponents-00");
  }
  
  m_sctmaterials = iAccessSvc->getRecordsetPtr("SCTMaterials",keySCT.tag(),keySCT.node());
  if(loadDefaults && m_sctmaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting SCTMaterials with default tag" <<endmsg;
    m_sctmaterials = iAccessSvc->getRecordsetPtr("SCTMaterials","SCTMaterials-00");
  }
  
  // --- TRT materials
  DecodeVersionKey keyTRT(iGeoModel, "TRT");
  m_trtmatcomponents = iAccessSvc->getRecordsetPtr("TrtMatComponents",keyTRT.tag(),keyTRT.node());
  if(loadDefaults && m_trtmatcomponents->size()==0)	{
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting TrtMatComponents with default tag" <<endmsg;
    m_trtmatcomponents = iAccessSvc->getRecordsetPtr("TrtMatComponents","TrtMatComponents-00");
  }
  m_trtmaterials = iAccessSvc->getRecordsetPtr("TrtMaterials",keyTRT.tag(),keyTRT.node());
  if(loadDefaults && m_trtmaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting TrtMaterials with default tag" <<endmsg;
    m_trtmaterials = iAccessSvc->getRecordsetPtr("TrtMaterials","TrtMaterials-00");
  }
  
  // --- InDet common materials
  DecodeVersionKey keyInDet(iGeoModel, "InnerDetector");
  m_indetmatcomponents = iAccessSvc->getRecordsetPtr("InDetMatComponents",keyInDet.tag(),keyInDet.node());
  if(loadDefaults && m_indetmatcomponents->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting InDetMatComponents with default tag" <<endmsg;
    m_indetmatcomponents = iAccessSvc->getRecordsetPtr("InDetMatComponents","InDetMatComponents-00");
  }
  
  m_indetmaterials = iAccessSvc->getRecordsetPtr("InDetMaterials",keyInDet.tag(),keyInDet.node());
  if(loadDefaults && m_indetmaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting InDetMaterials with default tag" <<endmsg;
    m_indetmaterials = iAccessSvc->getRecordsetPtr("InDetMaterials","InDetMaterials-00");
  }
  
  // --- LAr materials
  DecodeVersionKey keyLAr(iGeoModel, "LAr");    
  m_larmatcomponents = iAccessSvc->getRecordsetPtr("LArMatComponents",keyLAr.tag(),keyLAr.node());
  if(loadDefaults && m_larmatcomponents->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting LArMatComponents with default tag" <<endmsg;
    m_larmatcomponents = iAccessSvc->getRecordsetPtr("LArMatComponents","LArMatComponents-00");
  }
  m_larmaterials = iAccessSvc->getRecordsetPtr("LArMaterials",keyLAr.tag(),keyLAr.node());
  if(loadDefaults && m_larmaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting LArMaterials with default tag" <<endmsg;
    m_larmaterials = iAccessSvc->getRecordsetPtr("LArMaterials","LArMaterials-00");
  }
  
  // --- Tile materials
  DecodeVersionKey keyTile(iGeoModel, "TileCal");    
  m_tilematcomponents = iAccessSvc->getRecordsetPtr("TileMatComponents",keyTile.tag(),keyTile.node());
  if(loadDefaults && m_tilematcomponents->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting TileMatComponents with default tag" <<endmsg;
    m_tilematcomponents = iAccessSvc->getRecordsetPtr("TileMatComponents","TileMatComponents-00");
  }
  m_tilematerials = iAccessSvc->getRecordsetPtr("TileMaterials",keyTile.tag(),keyTile.node());
  if(loadDefaults && m_tilematerials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting TileMaterials with default tag" <<endmsg;
    m_tilematerials = iAccessSvc->getRecordsetPtr("TileMaterials","TileMaterials-00");
  }
  
  // --- Muon
  DecodeVersionKey keyMuon(iGeoModel, "MuonSpectrometer");
  m_muomatcomponents = iAccessSvc->getRecordsetPtr("MUOMatComponents",keyMuon.tag(),keyMuon.node());
  if(loadDefaults && m_muomatcomponents->size()==0)	{
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting MUOMatComponents with default tag" <<endmsg;
    m_muomatcomponents = iAccessSvc->getRecordsetPtr("MUOMatComponents","MUOMatComponents-00");
  }
  m_muomaterials = iAccessSvc->getRecordsetPtr("MUOMaterials",keyMuon.tag(),keyMuon.node());
  if(loadDefaults && m_muomaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting MUOMaterials with default tag" <<endmsg;
    m_muomaterials = iAccessSvc->getRecordsetPtr("MUOMaterials","MUOMaterials-00");  
  }
  m_shieldmatcomponents = iAccessSvc->getRecordsetPtr("ShieldMatComponents",keyMuon.tag(),keyMuon.node());
  if(loadDefaults && m_shieldmatcomponents->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting ShieldMatComponents with default tag" <<endmsg;
    m_shieldmatcomponents = iAccessSvc->getRecordsetPtr("ShieldMatComponents","ShieldMatComponents-00");
  }
  m_shieldmaterials = iAccessSvc->getRecordsetPtr("ShieldMaterials",keyMuon.tag(),keyMuon.node());
  if(loadDefaults && m_shieldmaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting ShieldMaterials with default tag" <<endmsg;
    m_shieldmaterials = iAccessSvc->getRecordsetPtr("ShieldMaterials","ShieldMaterials-00");
  }
  m_toromatcomponents = iAccessSvc->getRecordsetPtr("ToroMatComponents",keyMuon.tag(),keyMuon.node());
  if(loadDefaults && m_toromatcomponents->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting ToroMatComponents with default tag" <<endmsg;
    m_toromatcomponents =	iAccessSvc->getRecordsetPtr("ToroMatComponents","ToroMatComponents-00");
  }
  m_toromaterials = iAccessSvc->getRecordsetPtr("ToroMaterials",keyMuon.tag(),keyMuon.node());
  if(loadDefaults && m_toromaterials->size()==0) {
    if(log.level()<=MSG::WARNING)
      log << MSG::WARNING << " Getting ToroMaterials with default tag" <<endmsg; 
    m_toromaterials = iAccessSvc->getRecordsetPtr("ToroMaterials","ToroMaterials-00");
  }
  return StatusCode::SUCCESS;
}

// Destructor:
RDBMaterialManager::~RDBMaterialManager() {
	
  // Unreference the materials:
  for (auto &p : m_materialMap) {
    p.second->unref();
  }
	 	
  // Unreference the elements:
  for (GeoElement *elt : m_elementVector) {
    elt->unref();
  }
}

GeoMaterial* RDBMaterialManager::searchMaterialMap(const std::string & name) const
{
  MsgStream log(Athena::getMessageSvc(), "GeoModelSvc::RDBMaterialManager"); 
	
  std::map< std::string, GeoMaterial * >::const_iterator m   = m_materialMap.find(std::string(name));
  std::map< std::string, GeoMaterial * >::const_iterator end = m_materialMap.end();
  if (m!=end) {
    if(log.level()==MSG::VERBOSE)
      log << MSG::VERBOSE << " ***** in searchMaterialMap(): search sucess "  << endmsg;	
    return (*m).second;
  }
  else {
    if(log.level()==MSG::VERBOSE)
      log << MSG::VERBOSE << " ***** in searchMaterialMap(): search fail "  << endmsg;	
    return NULL;
  }
}


GeoElement *RDBMaterialManager::searchElementVector(const std::string & name)  const
{ 
  MsgStream log(Athena::getMessageSvc(), "GeoModelSvc::RDBMaterialManager"); 
	
  NameEquals matchByName(name);
  std::vector<GeoElement *>::const_iterator e=std::find_if(m_elementVector.begin(), m_elementVector.end(),matchByName);
  	
  if (e!=m_elementVector.end()) {	
    if(log.level()==MSG::VERBOSE)    		
      log << MSG::VERBOSE << " ***** in searchElementVector() search succes "  << endmsg;	
    return *e;
  }
  else {
    if(log.level()==MSG::VERBOSE)
      log << MSG::VERBOSE << " ***** in searchElementVector() search fail "  << endmsg;	
    return NULL;
  }
}


GeoElement *RDBMaterialManager::searchElementVector(const unsigned int atomicNumber) const
{ 
  MsgStream log(Athena::getMessageSvc(), "GeoModelSvc::RDBMaterialManager"); 
	
  NumberEquals matchByNumber(atomicNumber);
  std::vector<GeoElement *>::const_iterator e=std::find_if(m_elementVector.begin(), m_elementVector.end(), matchByNumber);
  	
  if (e!=m_elementVector.end()) {
    if(log.level()==MSG::VERBOSE)  		
      log << MSG::VERBOSE << " ***** in searchElementVector(atomicNumber) search succes "  << endmsg;
    return *e;
  }
  else {
    if(log.level()==MSG::VERBOSE)
      log << MSG::VERBOSE << " ***** in searchElementVector(atomicNumber) search succes "  << endmsg;
    return NULL;
  }
}

const GeoMaterial*  RDBMaterialManager:: getMaterial(const std::string &name) {

  unsigned int  ind, com_ind;
	
  std::string material_name;
  std::string tmp_name;
  long 	    material_id = 0;
  double      material_density = 0;
	
	
  std::string component_name;
  double      component_fraction;
  int 	    component_id;
		
  std::string detector;
  std::string tmp_det;
  std::string data_id;
	
	
  std::string matcomponents_table;

  [[maybe_unused]] static const bool specialMaterialsDone = [this]() {
    buildSpecialMaterials();
    return true;
  }();

  GeoMaterial* pmaterial;

  const GeoElement*  p_com_element;
	
  IRDBRecordset_ptr tmp_materials;
  IRDBRecordset_ptr tmp_matcomponents;
	
  MsgStream log(Athena::getMessageSvc(), "GeoModelSvc::RDBMaterialManager"); 
  if(log.level()<=MSG::DEBUG) 
    log << MSG::DEBUG  << " ***** getMaterial( ): "  << name << endmsg;	

  pmaterial = NULL;
  pmaterial = searchMaterialMap( name);
  if (pmaterial!= NULL) 
      return pmaterial;

  if(boost::starts_with(name, "std"))
    {
      detector = "std";
      tmp_materials = m_stdmaterials;
      tmp_matcomponents = m_stdmatcomponents;
      data_id = "STDMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "trt"))
    {
      detector = "trt";
      tmp_materials = m_trtmaterials;
      tmp_matcomponents = m_trtmatcomponents;
      data_id = "TRTMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "LAr"))
    {
      detector = "LAr";
      tmp_materials = m_larmaterials;
      tmp_matcomponents = m_larmatcomponents;
      data_id = "LARMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "muo"))
    {
      detector = "muo";
      tmp_materials = m_muomaterials;
      tmp_matcomponents = m_muomatcomponents;
      data_id = "MUOMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "pixtb"))
    {
      detector = "pixtb";
      tmp_materials = m_pixtbmaterials;
      tmp_matcomponents = m_pixtbmatcomponents;
      data_id = "PIXELTBMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "pix"))
    {
      detector = "pix";
      tmp_materials = m_pixmaterials;
      tmp_matcomponents = m_pixmatcomponents;
      data_id = "PIXMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "sct"))
    {
      detector = "sct";
      tmp_materials = m_sctmaterials;
      tmp_matcomponents = m_sctmatcomponents;
      data_id = "SCTMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "indet"))
    {
      detector = "indet";
      tmp_materials = m_indetmaterials;
      tmp_matcomponents = m_indetmatcomponents;
      data_id = "INDETMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "shield"))
    {
      detector = "shield";
      tmp_materials = m_shieldmaterials;
      tmp_matcomponents = m_shieldmatcomponents;
      data_id = "SHIELDMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "tile"))
    {
      detector = "tile";
      tmp_materials = m_tilematerials;
      tmp_matcomponents = m_tilematcomponents;
      data_id = "TILEMATERIALS_DATA_ID";
    }
  else if(boost::starts_with(name, "toro"))
    {
      detector = "toro";
      tmp_materials = m_toromaterials;
      tmp_matcomponents = m_toromatcomponents;
      data_id = "TOROMATERIALS_DATA_ID";
    }
  else {return 0 ;}
	
  for( ind = 0; ind < tmp_materials->size(); ind++)
    {
      const IRDBRecord* rec = (*tmp_materials)[ind];
      tmp_name = detector+"::"+rec->getString("NAME");
		
      if( name == tmp_name){
	material_name  =detector+"::"+rec->getString("NAME");
	material_id = rec->getLong(data_id);
	material_density = rec->getDouble("DENSITY");
        		
	if(log.level()<=MSG::DEBUG)
	  log << MSG::DEBUG  << " ***** Material: name id density: "  << material_name <<" " << material_id <<" "<< material_density << endmsg;	
	break;
      }
    }
		
  if (ind == tmp_materials->size()) 
      return NULL;

  pmaterial = new GeoMaterial( material_name,material_density * (GeoModelKernelUnits::gram / Gaudi::Units::cm3));

  bool firstComponent = true;
  bool hasSubMaterial = false;
  bool calculateFraction = false;
  double totalFraction = 0.;
  std::vector <const GeoElement*> elementComponents;
  std::vector <double>        elementFractions;

  for(  com_ind = 0; com_ind <tmp_matcomponents->size(); com_ind++)
    {
      const IRDBRecord* com_rec = (*tmp_matcomponents)[com_ind];
		
      component_id = com_rec->getLong("MATERIAL_ID");
      if( component_id == material_id)
	{
	  component_name = com_rec->getString("COMPNAME");
	  component_fraction = com_rec->getDouble("FRACTION");
			
	  if(firstComponent)
	  {
	    firstComponent = false;
	    if(component_fraction>=1.)
	      calculateFraction = true;
	  }

	  if( CheckElement( component_name) == 1)
	    {
	      p_com_element = getElement(component_name);

	      if(calculateFraction)
	      {
		totalFraction += component_fraction*p_com_element->getA();
		elementComponents.push_back(p_com_element);
		elementFractions.push_back(component_fraction);
	      }
	      else
		pmaterial->add( p_com_element, component_fraction);
										
	    }
	  else{
	    hasSubMaterial = true;
	    const GeoMaterial* p_com_material = getMaterial(component_name);
	    pmaterial->add(p_com_material, component_fraction);
			
	  }		
	}
    }    

  if(calculateFraction && hasSubMaterial && elementComponents.size()>0)
    std::cerr << material_name << " description should be changed. Please indicate the exact fraction for elements\n";

  if(calculateFraction && !elementComponents.empty()) {
    double inv_totalFraction = 1. / totalFraction;
    for(unsigned i=0; i<elementComponents.size(); i++)
      pmaterial->add(elementComponents[i],elementFractions[i]*elementComponents[i]->getA() * inv_totalFraction);
  }

  // a table to keep the memory allocation, and easy for delete
  addMaterial(detector,pmaterial);
	
  return pmaterial;
}


const GeoElement *RDBMaterialManager::getElement(const std::string & name) {
	
  unsigned int ind;

  std::string element_name;
  std::string element_symbol;
  std::string tmp_name;
	
  double      element_a;
  double      element_z;
	
  GeoElement *pelement;

  pelement = NULL;
  pelement = searchElementVector( name);
  if (pelement != NULL) 
      return pelement;

  MsgStream log(Athena::getMessageSvc(), "GeoModelSvc::RDBMaterialManager"); 
  if(log.level()==MSG::VERBOSE)
    log << MSG::VERBOSE << " ***** getElement(): " << name  <<endmsg;

  for(ind = 0; ind < m_elements->size(); ind++)
    {
      const IRDBRecord* rec = (*m_elements)[ind];
		
      tmp_name = rec->getString("NAME");
	
      if( name == tmp_name)
	{ 
	  element_name   = rec->getString("NAME");
	  element_symbol = rec->getString("SYMBOL");
	  element_a = rec->getDouble("A");
	  element_z = rec->getDouble("Z");
                	
	  pelement = new GeoElement( element_name , element_symbol  ,element_z , element_a *(GeoModelKernelUnits::gram/Gaudi::Units::mole));

	  // a table to keep the memory allocation, and easy for delete 
	  pelement->ref();
	  m_elementVector.push_back( pelement);
			
	  break;
	}
    }
  if (ind == m_elements->size()) 		return NULL;
	
  return pelement;

}


const GeoElement *RDBMaterialManager::getElement(unsigned int atomicNumber) {

  unsigned int ind;

  std::string element_name;
  std::string element_symbol;
	
  double      element_a;
  double      element_z;
	
  GeoElement* pelement(0);

  MsgStream log(Athena::getMessageSvc(), "GeoModelSvc::RDBMaterialManager"); 
  if(log.level()==MSG::VERBOSE)
    log << MSG::VERBOSE << " ***** const getElement(atomicNumber) const : " << atomicNumber <<endmsg;	

  for(ind = 0; ind < m_elements->size(); ind++)
    {
      const IRDBRecord* rec = (*m_elements)[ind];
		
      if(atomicNumber == rec->getDouble("A"))
	{ 
	  element_name   = rec->getString("NAME");
	  element_symbol = rec->getString("SYMBOL");
	  element_a = rec->getDouble("A");
	  element_z = rec->getDouble("Z");
                	
	  pelement = new GeoElement( element_name , element_symbol  ,element_z , element_a *(GeoModelKernelUnits::gram/Gaudi::Units::mole));

	  // a table to keep the memory allocation, and easy for delete 
	  pelement->ref();
	  m_elementVector.push_back( pelement);
			
	  break;
	}
    }
  if (ind == m_elements->size()) 	return NULL;
	
  return pelement;
}

void RDBMaterialManager::addMaterial(const std::string & /*space*/, GeoMaterial *material) {
	
  MsgStream log(Athena::getMessageSvc(), "GeoModelSvc::RDBMaterialManager"); 
  if(log.level()==MSG::VERBOSE)
    log << MSG::VERBOSE << " ***** RDBMaterialManager::addMaterial() "<<endmsg;
	
  std::string key = std::string(material->getName());
  // Check whether we already have materials with the same space::name defined
  if(m_materialMap.find(key)!=m_materialMap.end())
    log << MSG::WARNING << " Attempt to redefine material " << key << "!. The existing instance is kept. Please choose another name for new material" << endmsg;
  else {
    material->lock();             
    material->ref();
    m_materialMap[key]=material;
  }
}

StoredMaterialManager::MaterialMapIterator RDBMaterialManager::begin() const
{
  return m_materialMap.begin();
}

StoredMaterialManager::MaterialMapIterator RDBMaterialManager::end() const
{
  return m_materialMap.end();
}

size_t RDBMaterialManager::size()
{
  return m_materialMap.size();
}

std::ostream &  RDBMaterialManager::printAll(std::ostream & o) const 
{
  o << "============Material Manager Element List========================" << std::endl;
  for (GeoElement* elt : m_elementVector)
    {
      o << elt->getSymbol() << '\t' << elt->getZ() <<  '\t' << elt->getA() * (Gaudi::Units::mole / GeoModelKernelUnits::gram) << '\t' << elt->getName() << std::endl;
    }

  for (const auto& p : m_materialMap)
    {
      o << "Material: " << p.first <<  " Density " << p.second->getDensity() * (Gaudi::Units::cm3 / GeoModelKernelUnits::gram)  << std::endl;
      for (size_t i = 0; i< p.second->getNumElements();i++) 
	{
	  o <<" ***** ***** "<< int (p.second->getFraction(i)*100) << "% \t"  << p.second->getElement(i)->getName() << std::endl;
	}
    }
  	  	
  return o;
}

void RDBMaterialManager::buildSpecialMaterials()
{
  // Create special materials
  GeoElement* ethElement = new GeoElement("Ether","ET",500.0,0.0);
  ethElement->ref();
  m_elementVector.push_back(ethElement);  
  GeoMaterial* ether = new GeoMaterial("special::Ether",0.0);	
  ether->add(ethElement,1.);
  addMaterial("special",ether);
  // "Alternative" assembly material
  GeoMaterial* hu = new GeoMaterial("special::HyperUranium",0.0);	
  hu->add(ethElement,1.);
  addMaterial("special",hu);
}
