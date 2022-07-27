/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack



//
// includes
//

#include <EgammaAnalysisAlgorithms/EgammaIsolationCorrectionAlg.h>

//
// method implementations
//

namespace CP
{
  EgammaIsolationCorrectionAlg ::
  EgammaIsolationCorrectionAlg (const std::string& name, 
                     ISvcLocator* pSvcLocator)
    : AnaAlgorithm (name, pSvcLocator)
    , m_isolationCorrectionTool ("CP::IsolationCorrectionTool", this)
  {
    declareProperty ("isolationCorrectionTool", m_isolationCorrectionTool, "the smearing tool we apply");
  }



  StatusCode EgammaIsolationCorrectionAlg ::
  initialize ()
  {
    ANA_CHECK (m_isolationCorrectionTool.retrieve());
    ANA_CHECK (m_egammaHandle.initialize (m_systematicsList));
    ANA_CHECK (m_preselection.initialize (m_systematicsList, m_egammaHandle, SG::AllowEmpty));
    ANA_CHECK (m_systematicsList.addSystematics (*m_isolationCorrectionTool));
    ANA_CHECK (m_systematicsList.initialize());
    ANA_CHECK (m_outOfValidity.initialize());
    return StatusCode::SUCCESS;
  }



  StatusCode EgammaIsolationCorrectionAlg ::
  execute ()
  {
    for (const auto& sys : m_systematicsList.systematicsVector())
    {
      ANA_CHECK (m_isolationCorrectionTool->applySystematicVariation (sys));
      xAOD::EgammaContainer *egammas = nullptr;
      ANA_CHECK (m_egammaHandle.getCopy (egammas, sys));
      for (xAOD::Egamma *egamma : *egammas)
      {
        if (m_preselection.getBool (*egamma, sys))
        {
          ANA_CHECK_CORRECTION (m_outOfValidity, *egamma, m_isolationCorrectionTool->applyCorrection (*egamma));
        }
      }
    }
    return StatusCode::SUCCESS;
  }
}
