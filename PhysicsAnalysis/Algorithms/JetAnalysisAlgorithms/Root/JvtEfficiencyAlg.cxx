/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack


//
// includes
//

#include <JetAnalysisAlgorithms/JvtEfficiencyAlg.h>

#include <SelectionHelpers/SelectionHelpers.h>

//
// method implementations
//

namespace CP
{
  JvtEfficiencyAlg ::
  JvtEfficiencyAlg (const std::string& name, 
                     ISvcLocator* pSvcLocator)
    : AnaAlgorithm (name, pSvcLocator)
    , m_efficiencyTool ("", this)
    , m_truthJetsName("AntiKt4TruthJets")
  {
    declareProperty ("efficiencyTool", m_efficiencyTool, "the efficiency tool we apply");
    declareProperty ("dofJVT", m_dofJVT, "differenciate between JVT and fJVT");
    declareProperty ("fJVTStatus", m_fJVTStatus, "the decoration for the fJVT status");
    declareProperty ("skipBadEfficiency", m_skipBadEfficiency, "whether to skip efficiency calculation if the selection failed");
    declareProperty ("truthJetCollection", m_truthJetsName, "the truth jet collection to use for truth tagging");
  }



  StatusCode JvtEfficiencyAlg ::
  initialize ()
  {
    if (m_dofJVT && m_fJVTStatus.empty())
    {
      ANA_MSG_ERROR ("fJVTStatus decoration needs to be configured when running fJVT");
      return StatusCode::FAILURE;
    }

    ANA_CHECK (m_efficiencyTool.retrieve());
    ANA_CHECK (m_jetHandle.initialize (m_systematicsList));
    ANA_CHECK (m_preselection.initialize (m_systematicsList, m_jetHandle, SG::AllowEmpty));
    ANA_CHECK (m_selectionHandle.initialize (m_systematicsList, m_jetHandle, SG::AllowEmpty));
    ANA_CHECK (m_scaleFactorDecoration.initialize (m_systematicsList, m_jetHandle, SG::AllowEmpty));
    ANA_CHECK (m_systematicsList.addSystematics (*m_efficiencyTool));
    ANA_CHECK (m_systematicsList.initialize());
    ANA_CHECK (m_outOfValidity.initialize());

    if (m_dofJVT && !m_fJVTStatus.empty())
      ANA_CHECK (makeSelectionReadAccessor (m_fJVTStatus, m_fJVTStatusAccessor));

    return StatusCode::SUCCESS;
  }



  StatusCode JvtEfficiencyAlg ::
  execute ()
  {
    for (const auto& sys : m_systematicsList.systematicsVector())
    {
      ANA_CHECK (m_efficiencyTool->applySystematicVariation (sys));
      const xAOD::JetContainer *jets = nullptr;
      ANA_CHECK (m_jetHandle.retrieve (jets, sys));

      const xAOD::JetContainer *truthjets = nullptr;
      if(!m_truthJetsName.empty()) {
        ANA_CHECK(evtStore()->retrieve(truthjets,m_truthJetsName));
        ANA_CHECK(m_efficiencyTool->tagTruth(jets,truthjets));
      }

      for (const xAOD::Jet *jet : *jets)
      {
        if (m_preselection.getBool (*jet, sys))
        {
          bool goodJet = true;
          if (m_selectionHandle || m_skipBadEfficiency)
          {
            goodJet = m_dofJVT ? m_fJVTStatusAccessor->getBool (*jet) : m_efficiencyTool->passesJvtCut (*jet);
            if (m_selectionHandle)
              m_selectionHandle.setBool (*jet, goodJet, sys);
          }
          if (m_scaleFactorDecoration)
          {
            float sf = 1;
            if (goodJet) {
              ANA_CHECK_CORRECTION (m_outOfValidity, *jet, m_efficiencyTool->getEfficiencyScaleFactor (*jet, sf));
            } else if (!m_skipBadEfficiency) {
              ANA_CHECK_CORRECTION (m_outOfValidity, *jet, m_efficiencyTool->getInefficiencyScaleFactor (*jet, sf));
            }
            m_scaleFactorDecoration.set (*jet, sf, sys);
          }
        } else {
          if (m_selectionHandle)
            m_selectionHandle.setBool (*jet, false, sys);

          if (m_scaleFactorDecoration)
            m_scaleFactorDecoration.set (*jet, invalidScaleFactor(), sys);
        }
      }
    }

    return StatusCode::SUCCESS;
  }
}
