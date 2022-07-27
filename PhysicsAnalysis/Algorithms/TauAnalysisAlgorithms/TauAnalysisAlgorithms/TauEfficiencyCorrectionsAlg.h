/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack


#ifndef TAU_ANALYSIS_ALGORITHMS__TAU_EFFICIENCY_CORRECTIONS_ALG_H
#define TAU_ANALYSIS_ALGORITHMS__TAU_EFFICIENCY_CORRECTIONS_ALG_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <TauAnalysisTools/ITauEfficiencyCorrectionsTool.h>
#include <SelectionHelpers/OutOfValidityHelper.h>
#include <SelectionHelpers/SysReadSelectionHandle.h>
#include <SystematicsHandles/SysReadHandle.h>
#include <SystematicsHandles/SysWriteDecorHandle.h>
#include <SystematicsHandles/SysListHandle.h>
#include <SystematicsHandles/SysReadHandle.h>
#include <xAODTau/TauJetContainer.h>

namespace CP
{
  /// \brief an algorithm for calling \ref ITauEfficiencyCorrectionsTool

  class TauEfficiencyCorrectionsAlg final : public EL::AnaAlgorithm
  {
    /// \brief the standard constructor
  public:
    TauEfficiencyCorrectionsAlg (const std::string& name, 
                                   ISvcLocator* pSvcLocator);


  public:
    StatusCode initialize () override;

  public:
    StatusCode execute () override;
    


    /// \brief the smearing tool
  private:
    ToolHandle<TauAnalysisTools::ITauEfficiencyCorrectionsTool> m_efficiencyCorrectionsTool;

    /// \brief the systematics list we run
  private:
    SysListHandle m_systematicsList {this};

    /// \brief the tau collection we run on
  private:
    SysReadHandle<xAOD::TauJetContainer> m_tauHandle {
      this, "taus", "TauJets", "the tau collection to run on"};

    /// \brief the preselection we apply to our input
  private:
    SysReadSelectionHandle m_preselection {
      this, "preselection", "", "the preselection to apply"};

    /// \brief the helper for OutOfValidity results
  private:
    OutOfValidityHelper m_outOfValidity {this};

    /// \brief the decoration for the tau scale factor
  private:
    SysWriteDecorHandle<float> m_scaleFactorDecoration {
      this, "scaleFactorDecoration", "", "the decoration for the tau efficiency scale factor"};
  };
}

#endif
