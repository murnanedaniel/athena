/*
Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/ 
#ifndef InDet_PixelServicesTool_H
#define InDet_PixelServicesTool_H

#include "PixelInterfaces/IPixelServicesTool.h"
#include "PixelGeoModel/PixelGeoBuilder.h"
#include "PixelServicesTool/ServiceDynamicBuilder.h"

#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaKernel/IOVSvcDefs.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ServiceHandle.h"

//#include "InDetGeoModelStdUtils/GeoEnvelopeStd.h"

class InDetMaterialManager;
class OraclePixGeoAccessor;
class MsgStreamMember;

class IInDetServMatBuilderTool;

//class GeoPixelServices;
class PixelGeoBuilderBasics;

//namespace InDet {  

class ServiceStaticBuilder;
class ServiceDynamicBuilder;
class ServiceDynamicRouteBuilder;
class GeoServiceAssembly;
class GeoSimpleObject;
namespace InDetDD{
  class TubeZone;
}

class PixelServicesTool : virtual public IPixelServicesTool, public AthAlgTool {
  
 public:
  
  PixelServicesTool(const std::string&,const std::string&,const IInterface*);
  
  /** default destructor */
  virtual ~PixelServicesTool();
  
  virtual StatusCode initialize();
  /*    virtual StatusCode create(); */
  virtual StatusCode finalize();
  
  /*     // Register callback function on ConDB object */
  /*    virtual StatusCode registerCallback( StoreGateSvc* detStore ); */
     
  //    void build( const OraclePixGeoAccessor* accessor, InDetMaterialManager* matMgr );
     void buildServices(const PixelGeoBuilderBasics*, std::vector<InDetDD::TubeZone*> v = std::vector<InDetDD::TubeZone*>() );
		 
  void buildAndPlace(const std::string & region, GeoPhysVol * parent, double zcenter = 0, 
		     std::vector<std::string> svcList = std::vector<std::string>(),
		     bool bStatic = true, bool bDynamic = true);
  void buildAndPlace(const std::string & region, GeoFullPhysVol * parent, double zcenter = 0, 
		     std::vector<std::string> svcList = std::vector<std::string>(),
		     bool bStatic = true, bool bDynamic = true);
  double computeRmin(const std::string & region, std::vector<std::string> svcList = std::vector<std::string>()) const;
  double computeRmax(const std::string & region, std::vector<std::string> svcList = std::vector<std::string>()) const;
  double computeZmin(const std::string & region, std::vector<std::string> svcList = std::vector<std::string>()) const;
  double computeZmax(const std::string & region, std::vector<std::string> svcList = std::vector<std::string>()) const;
  std::string getLayerModuleMaterialName(int, int) const;
  std::string getLayerModuleMaterialName(int, std::vector<int>) const;
  std::string getLayerStaveModuleMaterialName(int, int, int) const;
  std::string getLayerStaveModuleMaterialName(int, int, std::vector<int>) const;

  void resetServices();
  bool svcRouteAuto() const { return (m_dynServices!=0&&!m_bSvcDynAutomated); }
  //    std::vector<InDet::GeoServiceAssembly* > getServiceAssemblies();
  //    std::vector<InDet::GeoSimpleObject* > getServiceObjects();

 InDetMaterialManager* materialMgr() const {return m_matMgr;}

 double getMaterialFudgeModuleSvc(int iLayer) const;
 double getMaterialFudgeSvcEc(int iLayer) const;
 ServiceDynamicBuilder::SvcEcMaterialFudges getMaterialFudgesSvcEc() const;
  
 private:
  
  ServiceStaticBuilder * m_pixServices;
  mutable ServiceDynamicBuilder * m_dynServices;
  ServiceDynamicBuilder::SvcEcMaterialFudges m_ECmatrialFudges;
  //  mutable ServiceDynamicRouteBuilder * m_dynRouteServices;
  //  ToolHandle< IInDetServMatBuilderTool > m_serviceBuilderTool;

  bool m_bReadSvcFromDB;
  bool m_bSetupBarrelModuleMaterial;
  bool m_bSvcDynAutomated;
  std::string m_svcListFromDB;
  
  mutable InDetMaterialManager* m_matMgr;
  mutable Athena::MsgStreamMember m_msg;
  const PixelGeoBuilderBasics* m_basics;

};

//}

#endif
