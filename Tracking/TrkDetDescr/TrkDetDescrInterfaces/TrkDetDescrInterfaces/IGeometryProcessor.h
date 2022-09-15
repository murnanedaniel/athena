/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// IGeometryProcessor.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKDETDESCRINTERFACES_IGEOMETRYPROCESSOR_H
#define TRKDETDESCRINTERFACES_IGEOMETRYPROCESSOR_H

// Gaudi
#include "GaudiKernel/IAlgTool.h"
//STL
#include <string>

#include "CxxUtils/checker_macros.h"
namespace Trk {

  class TrackingGeometry;
  class TrackingVolume;
  class Layer;
  class Surface;

  /** Interface ID for IGeometryProcessors*/  
  static const InterfaceID IID_IGeometryProcessor("IGeometryProcessor", 1, 0);
  
  /** @class IGeometryProcessor
  
       Interface class IGeometryProcessors
      @author Andreas.Salzburger@cern.ch
      @author Christos Anastopoulos MT fixes
    */
  class IGeometryProcessor : virtual public IAlgTool {
    
    public:
      /**Virtual destructor*/
      virtual ~IGeometryProcessor(){}
      
      /** AlgTool and IAlgTool interface methods */
      static const InterfaceID& interfaceID() { return IID_IGeometryProcessor; }

      /** Processor Action to work on TrackingGeometry& tgeo */
      virtual StatusCode process  (TrackingGeometry& tvol) const = 0;

      /** Processor Action to work on TrackingVolumes - the level is for the hierachy tree*/
      virtual StatusCode process  (TrackingVolume& tvol, size_t level=0) const = 0;
     
      /** Processor Action to work on Layers */
      virtual StatusCode process  (Layer& lay, size_t level=0) const = 0;

      /** Processor Action to work on Surfaces */
      virtual StatusCode process(Surface& surf, size_t level=0) const = 0;

  };

} // end of namespace

#endif // TRKDETDESCRINTERFACES_IGEOMETRYPROCESSOR_H
