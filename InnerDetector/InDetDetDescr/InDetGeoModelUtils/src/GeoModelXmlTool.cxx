/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "InDetGeoModelUtils/GeoModelXmlTool.h"

#include <GeoModelKernel/GeoPhysVol.h>
#include <GeoModelUtilities/DecodeVersionKey.h>
#include <GeoModelXml/Gmx2Geo.h>
#include <GeoModelXml/GmxInterface.h>
#include <PathResolver/PathResolver.h>
#include <RDBAccessSvc/IRDBRecord.h>
#include <RDBAccessSvc/IRDBRecordset.h>

#include <fstream>
#include <utility>

GeoModelXmlTool::GeoModelXmlTool(const std::string &type,
                                 const std::string &name,
                                 const IInterface *parent)
 : GeoModelTool(type, name, parent)
{
}

StatusCode GeoModelXmlTool::createBaseTool()
{
  ATH_CHECK(m_geoDbTagSvc.retrieve());
  ATH_CHECK(m_rdbAccessSvc.retrieve());

  return StatusCode::SUCCESS;
}

int GeoModelXmlTool::createTopVolume(GeoPhysVol* world, GmxInterface& gmxInterface, const std::string& vNode, const std::string& tableName) const
{
  int flags{};
  std::string gmxInput;

  if (m_gmxFilename.empty()) {
    ATH_MSG_INFO("Getting " << m_detectorName.value() << " GeoModelXml description from the geometry database");
    flags = 0x1; // Lowest bit ==> string; next bit implies gzip'd but we decided not to gzip
    //how to propagate these to here best...?
    gmxInput = getBlob(vNode,tableName);
    std::string dtdFile = '"' + PathResolver::find_file("GeoModelXml/geomodel.dtd", "DATAPATH") + '"';
    ATH_MSG_DEBUG("dtdFile = " << dtdFile);
    size_t index = gmxInput.find("\"geomodel.dtd\"");
    if (index != std::string::npos) {
      gmxInput.replace(index, 14, dtdFile);
    } else {
      throw std::runtime_error("GeoModeXmlDetectorFactoryBase::createTopVolume: Did not find string geomodel.dtd in the gmx input string.");
    }
  } else {
    flags = 0;
    gmxInput = PathResolver::find_file(m_gmxFilename, "DATAPATH");
    if (gmxInput.empty()) { // File not found
      std::string errMessage("GeoModeXmlDetectorFactoryBase::createTopVolume:: Unable to find file " + m_gmxFilename +
                             " with PathResolver; check filename and DATAPATH environment variable");
      throw std::runtime_error(errMessage);
    }
  }

  //optionally dump to local file for examination
  if(m_clobOutputFileName!="") {
      std::ofstream out(m_clobOutputFileName);
      out << gmxInput;
      out.close();
   }

  Gmx2Geo gmx2Geo(gmxInput, world, gmxInterface, flags);

  unsigned int nChildren = world->getNChildVols();

  bool foundVolume = false;
  int childIndex = -1;
  //find the appropriate volume in the hierarchy, to allow it to be set as the topVolume in
  //our detectorManager, by returning the index
  for (int iChild = nChildren - 1; iChild>=0; --iChild) {
    if (world->getNameOfChildVol(iChild) == m_detectorName) {
      // The * converts from a ConstPVLink to a reference to a GeoVPhysVol;
      // the & takes its address.
      foundVolume = true;
      childIndex = iChild;
      break;
    }
  }
  if(!foundVolume) ATH_MSG_ERROR("Couldn't find the volume in the world hierarchy!");
  return childIndex;
}

bool GeoModelXmlTool::isAvailable(const std::string& vNode, const std::string& tableName) const
{
  if (m_gmxFilename.empty()) {
    DecodeVersionKey versionKey(&*m_geoDbTagSvc, vNode);
    const std::string& versionTag  = versionKey.tag();
    const std::string& versionNode = versionKey.node();
    const std::string version = m_rdbAccessSvc->getChildTag(tableName, versionTag, versionNode);
    if (version.empty()) {
      return false;
    }
    ATH_MSG_INFO("Using " << version << " from " << versionNode << " tag " << versionTag);
  }

  return true;
}

std::string GeoModelXmlTool::getBlob(const std::string& vNode, const std::string& tableName) const
{
  DecodeVersionKey versionKey(&*m_geoDbTagSvc, vNode);
  const IRDBRecordset_ptr recordSet = m_rdbAccessSvc->getRecordsetPtr(tableName, versionKey.tag(), versionKey.node());
  if (!recordSet || recordSet->size() == 0) {
    ATH_MSG_FATAL("Unable to obtain " << vNode << " recordSet");
    throw std::runtime_error("Unable to obtain recordSet");
  }
  const IRDBRecord *record = (*recordSet)[0];
  std::string clobString = record->getString("XMLCLOB");
  return clobString;
}
