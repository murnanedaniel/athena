/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EGAMMAINTERFACES_IELECTRONSUPERCLUSTERBUILDER_H
#define EGAMMAINTERFACES_IELECTRONSUPERCLUSTERBUILDER_H

/** @class IELECTRONSUPERCLUSTERBULDER
  Interface for the tool Reconstruction/egamma/egammaTools/electronSuperClusterBuilder
  @author David Di Valentino david.di.valentino@cern.ch
CREATED  : 21/12/2013
MODIFIED : 
*/
// INCLUDE HEADER FILES:
#include "GaudiKernel/IAlgTool.h"

static const InterfaceID IID_IelectronSuperClusterBuilder("IelectronSuperClusterBuilder", 1, 0);

class IelectronSuperClusterBuilder : virtual public IAlgTool
{

 public:
  /** @brief Virtual destructor*/
  virtual ~IelectronSuperClusterBuilder() {};
  /** @brief AlgTool interface methods */
  static const InterfaceID& interfaceID();
  /** @brief initialize method*/
  virtual StatusCode initialize() = 0;
  virtual StatusCode finalize()   = 0;
  virtual StatusCode execute()    = 0; 
};

inline const InterfaceID& IelectronSuperClusterBuilder::interfaceID()
{
  return IID_IelectronSuperClusterBuilder;
}

#endif
