/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGINDETTOOLINTERFACES_ITRIGL2LAYERNUMBERTOOL_H
#define TRIGINDETTOOLINTERFACES_ITRIGL2LAYERNUMBERTOOL_H

#include "GaudiKernel/IAlgTool.h"
#include "TrigInDetPattRecoEvent/TrigInDetSiLayer.h"
#include <vector>


/** @class ITrigL2LayerNumberTool    
provides the abstract interface for the silicon layer number tool
@author D.Emeliyanov <http://consult.cern.ch/xwho>
*/

static const InterfaceID IID_ITrigL2LayerNumberTool("ITrigL2LayerNumberTool",1,0);

class ITrigL2LayerNumberTool : virtual public IAlgTool {
 public:
  /** other standard AlgTool methods */
  
  static const InterfaceID& interfaceID () {  //!< the Tool's interface
    return IID_ITrigL2LayerNumberTool; 
  }  	
  
  virtual int maxSiliconLayerNum() const = 0;
  virtual int offsetEndcapPixels() const = 0;
  virtual int offsetBarrelSCT() const    = 0;
  virtual int offsetEndcapSCT() const    = 0;
  virtual void report() const            = 0;//prints out the above

  virtual int maxNumberOfUniqueLayers() const  = 0;
  virtual const std::vector<short>* pixelLayers() const  = 0; // hashId addressable arrays of layer numbers
  virtual const std::vector<short>* sctLayers() const  = 0; // hashId addressable arrays of layer numbers
  virtual const std::vector<TrigInDetSiLayer>* layerGeometry() const = 0;
};

#endif
