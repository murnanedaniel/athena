/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef GEOMODELSVC_GEOMODELXMLTOOL_H
#define GEOMODELSVC_GEOMODELXMLTOOL_H

#include <GeoModelInterfaces/IGeoDbTagSvc.h>
#include <GeoModelUtilities/GeoModelTool.h>
#include <RDBAccessSvc/IRDBAccessSvc.h>

class GeoPhysVol;
class GeoVPhysVol;
class GmxInterface;

class GeoModelXmlTool : public GeoModelTool
{

public:
  GeoModelXmlTool(const std::string& type, const std::string& name, const IInterface* parent);

protected:
  bool isAvailable(const std::string& versionNode, const std::string& tableNode) const;
  std::string getBlob(const std::string& versionNode, const std::string& tableNode) const;
  // returns a position in the world volume hierarchy
  // -1 represents the volume named by m_detectorName not being found
  int createTopVolume(GeoPhysVol* worldVol, GmxInterface& interface, const std::string& versionNode, const std::string& tableNode) const;
  // version for detectors created within an envelope volume, returning a const GeoVPhysVol directly
  const GeoVPhysVol * createTopVolume(GeoPhysVol* worldVol, GmxInterface& interface, const std::string& versionNode, const std::string& tableNode, const std::string& containingDetector, const std::string& envelopeName) const;

  // method for derived classes to initialize the services needed here
  StatusCode createBaseTool();

  Gaudi::Property<std::string> m_gmxFilename{this, "GmxFilename", "", "The name of the local file to read the geometry from"};
  Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "ITkStrip", ""};
  ServiceHandle<IRDBAccessSvc> m_rdbAccessSvc{this, "RDBAccessSvc", "RDBAccessSvc", ""};
  ServiceHandle<IGeoDbTagSvc> m_geoDbTagSvc{this, "GeoDbTagSvc", "GeoDbTagSvc", ""};

private:

  Gaudi::Property<std::string> m_clobOutputFileName{this, "ClobOutputName", "", "Name of file to dump CLOB content to"};
  void createTopVolume_impl(GeoPhysVol* worldVol, GmxInterface& interface, const std::string& versionNode, const std::string& tableNode) const;

};

#endif // GEOMODELSVC_GEOMODELXMLTOOL_H
