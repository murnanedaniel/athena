/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef GENERATORFILTERS_XAODMISSINGETFILTER_H
#define GENERATORFILTERS_XAODMISSINGETFILTER_H

#include "GeneratorModules/GenFilter.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticleAuxContainer.h"


/// Filters on total missing energy from nus and LSPs
/// @author Seth Zenz, December 2005
class xAODMissingEtFilter:public GenFilter {
public:

  xAODMissingEtFilter(const std::string& name, ISvcLocator* pSvcLocator);
  virtual StatusCode filterEvent();

 private:

  double m_METmin;
  bool m_useHadronicNu;

};

#endif