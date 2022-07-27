/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef HGTD_GEOMODEL_HGTD_DETECTORFACTORY_H
#define HGTD_GEOMODEL_HGTD_DETECTORFACTORY_H
//
//    Main class for building the GeoModel geometry, and handle the GeometryManager and
//    DetectorManager.
//
#include "StoreGate/DataHandle.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "InDetGeoModelUtils/InDetDetectorFactoryBase.h"
#include "HGTD_ReadoutGeometry/HGTD_DetectorManager.h"

namespace InDetDD {
    class HGTD_ModuleDesign;
    class SiCommonItems;
}

class GeoTube ;
class GeoPhysVol;
class StoredMaterialManager;
class HGTD_GeoModelAthenaComps;

struct HgtdGeoParams {
  double rMid;
  double rOuter;
  double disk1Rotation;
  double disk2Rotation;
  double rowSpaceSide;
  double rowBacksideInnerShift;
  double rowBacksideOuterShift;
  double moduleSpaceInner;
  double moduleSpaceOuter;
  double flexSheetSpacing;
};

struct ModulePosition {
  double x;
  double y;
  double phiRotation;
  // below for backward compatibilty
  bool flipped;
  int row;
  int el_in_row;
};

struct GeoCylVolParams {
  std::string name;
  double rMin;
  double rMax;
  double zHalf;
  double zOffsetLocal;
  std::string material;
};

struct GeoBoxVolParams {
  std::string name;
  double xHalf;
  double yHalf;
  double zHalf;
  double zOffsetLocal;
  std::string material;
};

typedef std::array<std::vector<ModulePosition>,21> PositionsInQuadrant;


class HGTD_DetectorFactory : public InDetDD::DetectorFactoryBase {
public:
  HGTD_DetectorFactory(HGTD_GeoModelAthenaComps* athenaComps);
  virtual ~HGTD_DetectorFactory();

  // Creation of geometry:
  virtual void create(GeoPhysVol* world) override;

  // Access to the results:
  virtual const HGTD_DetectorManager* getDetectorManager() const override;

  void setPrintIdentifierDict( bool );

private:
  // Copy and assignments operations illegal and so are made private
  HGTD_DetectorFactory(HGTD_DetectorFactory &right);
  HGTD_DetectorFactory & operator=(HGTD_DetectorFactory &right);

  void initializeGeoParameters();
  GeoLogVol* buildEndcapLogicalVolume(bool isPositiveSide);
  GeoVPhysVol* build( const GeoLogVol* logicalEnvelope, bool bPos);

  InDetDD::HGTD_ModuleDesign* createHgtdDesign( double thickness );

  //  below 3 members prepare 3-ring vs 2-ring layout controlled implicitly by geomVersion
  std::array< PositionsInQuadrant, 4 > prepareLayersFromQuadrants( unsigned int ) ;
  PositionsInQuadrant prepareQuadrantsFromRows( int layer, unsigned int maxRow );
  std::string formModuleName( int layer, int quadrant, unsigned int maxrows, int row, int mod,
              ModulePosition module, double & myx, double & myy, double & myrot,
              int & phi, int & eta ) ;

  // 3-ring layout
  PositionsInQuadrant mirrorModulesInQuadrant( PositionsInQuadrant );
  std::vector< ModulePosition > prepareModulePositionsInRowThreeRing( int row, int back = 0 );
  int reorderRows( PositionsInQuadrant* quadrant );

  // 2-ring layout
  std::vector<ModulePosition> prepareModulePositionsInRowTwoRing(int row, bool back = false);

  void mirrorPositionsAroundYaxis(std::array< PositionsInQuadrant, 4 >& arr);

  HGTD_DetectorManager* m_detectorManager;
  HGTD_GeoModelAthenaComps* m_athComps;
  StoredMaterialManager* m_materialMgr;

  int m_geomVersion;

  // whether print number of modules per row for to the input for Identifier dictionary
  bool m_outputIdfr;

  std::map<std::string,GeoCylVolParams> m_cylVolPars;
  std::map<std::string,GeoBoxVolParams> m_boxVolPars;
  HgtdGeoParams m_hgtdPars;

  std::unique_ptr<const InDetDD::SiCommonItems> m_commonItems;
};

#endif // HGTD_GEOMODEL_HGTD_DETECTORFACTORY_H
