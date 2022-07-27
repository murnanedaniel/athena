/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

/// @author Miha Muskinja

//
// includes
//

#include <AsgAnalysisAlgorithms/EventSelectionByObjectFlagAlg.h>
#include <SystematicsHandles/SysFilterReporter.h>
#include <SystematicsHandles/SysFilterReporterCombiner.h>

//
// method implementations
//

namespace CP {

EventSelectionByObjectFlagAlg ::EventSelectionByObjectFlagAlg(
    const std::string &name, ISvcLocator *pSvcLocator)
    : AnaAlgorithm(name, pSvcLocator) {}

StatusCode EventSelectionByObjectFlagAlg ::initialize() {

    ANA_CHECK(m_particleHandle.initialize (m_systematicsList));
    ANA_CHECK(m_preselection.initialize (m_systematicsList, m_particleHandle, SG::AllowEmpty));
    ANA_CHECK(m_veto.initialize (m_systematicsList, m_particleHandle));
    ANA_CHECK(m_systematicsList.initialize());
    ANA_CHECK(m_filterParams.initialize());

    return StatusCode::SUCCESS;
}

StatusCode EventSelectionByObjectFlagAlg ::execute() {

    SysFilterReporterCombiner filterCombiner (m_filterParams, true);

    // loop over systematics
    for (const auto& sys : m_systematicsList.systematicsVector())
    {
        SysFilterReporter filter (filterCombiner, sys);

        // particle container
        const xAOD::IParticleContainer *particles = nullptr;
        ANA_CHECK(m_particleHandle.retrieve(particles, sys));

        // reject events with any particle passing the input selection
        for (const xAOD::IParticle *particle : *particles) {
           if (m_preselection.getBool(*particle, sys)) {
                if (m_veto.getBool (*particle, sys)) {
                    ATH_MSG_VERBOSE("Event failed.");
                    filter.setPassed (false);
                    break;
                }
            }
        }
    }

    return StatusCode::SUCCESS;
}

StatusCode EventSelectionByObjectFlagAlg ::finalize() {
    ANA_CHECK (m_filterParams.finalize());
    return StatusCode::SUCCESS;
}

} // namespace CP
