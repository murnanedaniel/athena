/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef GENERATORFILTERS_XAODTRUTHPARTICLESLIMMERPHOTON_H
#define GENERATORFILTERS_XAODTRUTHPARTICLESLIMMERPHOTON_H

#include "AthenaBaseComps/AthAlgorithm.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthMetaDataContainer.h"

/// @brief Algorithm to skim the xAOD truth particle container for xAOD photons filter
///
/// This algorithm is used to copy and skim the particles from the xAOD TruthParticles container,
/// keeping just relevant photons from the event.
/// The design of this class heavily mirrors the DerivationFramework::TruthCollectionMaker.
///
/// @author Jeff Dandoy <Jeff.Dandoy@cern.ch>
class xAODTruthParticleSlimmerPhoton : public AthAlgorithm
{
public:
    /// Regular algorithm constructor
    xAODTruthParticleSlimmerPhoton(const std::string &name, ISvcLocator *svcLoc);
    /// Function initialising the algorithm
    virtual StatusCode initialize();
    /// Function executing the algorithm
    virtual StatusCode execute();

private:
    /// The key for the output xAOD truth containers
    std::string m_xaodTruthParticleContainerNamePhoton;
    std::string m_xaodTruthParticleContainerName;
    std::string m_xaodTruthEventContainerName;
}; // class xAODTruthParticleSlimmerPhoton

#endif //GENERATORFILTERS_XAODTRUTHPARTICLESLIMMERPHOTON_H
