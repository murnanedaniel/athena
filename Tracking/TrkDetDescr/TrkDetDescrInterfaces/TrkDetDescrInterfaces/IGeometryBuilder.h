/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// IGeometryBuilder.hm (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKDETDESCRINTERFACES_IGEOMETRYBUILDER_H
#define TRKDETDESCRINTERFACES_IGEOMETRYBUILDER_H

// Gaudi
#include "GaudiKernel/IAlgTool.h"
// Trk - enum
#include "TrkDetDescrUtils/GeometrySignature.h"
#include "TrkSurfaces/Surface.h"
// STL
#include <vector>

#include "CxxUtils/checker_macros.h"

namespace Trk {

  class TrackingGeometry;
  class TrackingVolume;
  class Layer;

  /** @class IGeometryBuilder
    
    Interface class IGeometryBuilders,
    the GeometryBuilder inherits from this one.
    
    VolumeBounds can be given optionally to force a specific size/shape/boundary

    This interface class implements protected glue and surface
    exchange methods, that require friend rights to the classes
    
    @author Andreas.Salzburger@cern.ch
    */
  class IGeometryBuilder : virtual public IAlgTool {
    
    public:
    /// Creates the InterfaceID and interfaceID() method
    DeclareInterfaceID(IGeometryBuilder, 1, 0);

      /**Virtual destructor*/
      virtual ~IGeometryBuilder(){}

      /** TrackingGeometry Interface methode -
           - optionally a pointer to a tvol is given to wrap the TrackingGeometry around
           - optionally a vector of pointers to TrackingVolume is given for glueing
        */
      virtual std::unique_ptr<TrackingGeometry> trackingGeometry (TrackingVolume* tvol = 0) const = 0;

      /** The unique signature */
      virtual GeometrySignature geometrySignature() const = 0;
      
      protected:
      /** Protected method to register the Layer to the Surface */
      void associateLayer(const Layer& lay, Surface& sf) const
      {
        sf.associateLayer(lay);
      }
      
  };

} // end of namespace


#endif // TRKDETDESCRINTERFACES_IGEOMETRYBUILDER_H


