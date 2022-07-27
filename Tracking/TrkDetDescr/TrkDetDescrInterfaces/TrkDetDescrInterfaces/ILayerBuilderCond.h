/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// ILayerBuilderCond.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef TRKDETDESCRINTERFACES_ILAYERBUILDERCOND_H
#define TRKDETDESCRINTERFACES_ILAYERBUILDERCOND_H

// Gaudi
#include "GaudiKernel/EventContext.h"
#include "GaudiKernel/EventIDRange.h"
#include "GaudiKernel/IAlgTool.h"
//
#include "TrkSurfaces/Surface.h"
#include "StoreGate/WriteCondHandle.h"
// STL
#include <string>
#include <vector>

namespace Trk {

class CylinderLayer;
class DiscLayer;
class PlaneLayer;
class Layer;
class TrackingGeometry;

/** Interface ID for ILayerBuilderCond*/
static const InterfaceID IID_ILayerBuilderCond("ILayerBuilderCond", 1, 0);

/** @class ILayerBuilderCond
  Interface class ILayerBuilderConds
  It inherits from IAlgTool. The actual implementation of the AlgTool depends on
  the SubDetector, more detailed information should be found there.

  @author Andreas.Salzburger@cern.ch
  */
class ILayerBuilderCond : virtual public IAlgTool
{

public:
  /**Virtual destructor*/
  virtual ~ILayerBuilderCond() {}

  /** AlgTool and IAlgTool interface methods */
  static const InterfaceID& interfaceID() { return IID_ILayerBuilderCond; }

  /** LayerBuilder interface method - returning Barrel-like layers */
  virtual std::unique_ptr<const std::vector<CylinderLayer*> >
  cylindricalLayers(const EventContext& ctx,
                    SG::WriteCondHandle<TrackingGeometry>& whandle) const = 0;

  /** LayerBuilder interface method - returning Endcap-like layers */
  virtual std::unique_ptr<const std::vector<DiscLayer*> >
  discLayers(const EventContext& ctx,
             SG::WriteCondHandle<TrackingGeometry>& whandle) const = 0;

  /** LayerBuilder interface method - returning Planar-like layers */
  virtual std::unique_ptr<const std::vector<PlaneLayer*> >
  planarLayers(const EventContext& ctx,
               SG::WriteCondHandle<TrackingGeometry>& whandle) const = 0;

  /** Name identification */
  virtual const std::string& identification() const = 0;

  /** Validation Action:
      Can be implemented optionally, outside access to internal validation steps
   */
  virtual void validationAction() const {}

protected:
  /** Protected method to register the Layer to the Surface */
  void associateLayer(const Layer& lay, Surface& sf) const
  {
    sf.associateLayer(lay);
  }
};

} // end of namespace

#endif // TRKDETDESCRINTERFACES_ILAYERBUILDERCOND_H

