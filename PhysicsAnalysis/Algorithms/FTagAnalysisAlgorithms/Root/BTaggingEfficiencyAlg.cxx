/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack



//
// includes
//

#include <FTagAnalysisAlgorithms/BTaggingEfficiencyAlg.h>

//
// method implementations
//

namespace CP
{
  BTaggingEfficiencyAlg ::
  BTaggingEfficiencyAlg (const std::string& name, 
                     ISvcLocator* pSvcLocator)
    : AnaAlgorithm (name, pSvcLocator)
    , m_efficiencyTool ("BTaggingEfficiencyTool", this)
  {
    declareProperty ("efficiencyTool", m_efficiencyTool, "the calibration and smearing tool we apply");
    declareProperty ("onlyInefficiency", m_onlyInefficiency, "whether only to calculate inefficiencies");
  }



  StatusCode BTaggingEfficiencyAlg ::
  initialize ()
  {
    if (m_onlyInefficiency && m_selectionHandle)
    {
      ANA_MSG_ERROR ("can't specify both onlyInefficiency and selectionDecoration");
      return StatusCode::FAILURE;
    }

    if (m_scaleFactorDecoration.empty())
    {
      ANA_MSG_ERROR ("no scale factor decoration name set");
      return StatusCode::FAILURE;
    }

    ANA_CHECK (m_efficiencyTool.retrieve());
    ANA_CHECK (m_jetHandle.initialize (m_systematicsList));
    ANA_CHECK (m_preselection.initialize (m_systematicsList, m_jetHandle, SG::AllowEmpty));
    ANA_CHECK (m_selectionHandle.initialize (m_systematicsList, m_jetHandle, SG::AllowEmpty));
    ANA_CHECK (m_scaleFactorDecoration.initialize (m_systematicsList, m_jetHandle));
    ANA_CHECK (m_systematicsList.addSystematics (*m_efficiencyTool));
    ANA_CHECK (m_systematicsList.initialize());
    ANA_CHECK (m_outOfValidity.initialize());

    return StatusCode::SUCCESS;
  }



  StatusCode BTaggingEfficiencyAlg ::
  execute ()
  {
    for (const auto& sys : m_systematicsList.systematicsVector())
    {
      ANA_CHECK (m_efficiencyTool->applySystematicVariation (sys));
      const xAOD::JetContainer *jets = nullptr;
      ANA_CHECK (m_jetHandle.retrieve (jets, sys));
      for (const xAOD::Jet *jet : *jets)
      {
        if (m_preselection.getBool (*jet, sys))
        {
          float sf = 0;

          // The efficiency tool can calculate both efficiencies and
          // inefficiencies.  This setup can calculate either, or
          // both; in the case of the later a selection decoration is
          // used to decide whether to calculate efficiencies or
          // inefficiencies.
          //
          // Note that if you want to exclude jets from processing,
          // this selection accessor/decoration has nothing to do with
          // it.  You do the pre-selection via a view container like
          // for all the other CP algorithms.
          if (!m_onlyInefficiency && m_selectionHandle.getBool (*jet, sys))
          {
            ANA_CHECK_CORRECTION (m_outOfValidity, *jet, m_efficiencyTool->getScaleFactor (*jet, sf));
          } else
          {
            ANA_CHECK_CORRECTION (m_outOfValidity, *jet, m_efficiencyTool->getInefficiencyScaleFactor (*jet, sf));
          }
          m_scaleFactorDecoration.set (*jet, sf, sys);
        } else {
          m_scaleFactorDecoration.set (*jet, invalidScaleFactor(), sys);
        }
      }
    }
    return StatusCode::SUCCESS;
  }
}
