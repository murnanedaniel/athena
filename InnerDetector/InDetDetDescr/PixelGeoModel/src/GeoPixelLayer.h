/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef PIXELGEOMODEL_GEOPIXELLAYER_H
#define PIXELGEOMODEL_GEOPIXELLAYER_H

#include "GeoVPixelFactory.h"

class GeoPhysVol;
class GeoTransform;

class GeoPixelLayer : public GeoVPixelFactory {

 public:
  GeoPixelLayer(InDetDD::PixelDetectorManager* ddmgr,
                PixelGeometryManager* mgr);
  virtual GeoVPhysVol* Build() override;

  GeoPhysVol* getSupportA(){ return m_supportPhysA; }
  GeoPhysVol* getSupportC(){ return m_supportPhysC; }
  GeoVPhysVol* getSupportMidRing(){ return m_supportMidRing; }

  GeoTransform* getSupportTrfA(){ return m_xformSupportA; }
  GeoTransform* getSupportTrfC(){ return m_xformSupportC; }
  GeoTransform* getSupportTrfMidRing(){ return m_xformSupportMidRing; }

 private:
  GeoPhysVol *m_supportPhysA;
  GeoPhysVol *m_supportPhysC;
  GeoVPhysVol *m_supportMidRing;

  GeoTransform *m_xformSupportA;
  GeoTransform *m_xformSupportC;
  GeoTransform *m_xformSupportMidRing;

};

#endif // not PIXELGEOMODEL_GEOPIXELLAYER_H
