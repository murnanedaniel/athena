/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack


#ifndef F_TAG_ANALYSIS_ALGORITHMS__B_TAGGING_EFFICIENCY_ALG_H
#define F_TAG_ANALYSIS_ALGORITHMS__B_TAGGING_EFFICIENCY_ALG_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <FTagAnalysisInterfaces/IBTaggingEfficiencyTool.h>
#include <SelectionHelpers/OutOfValidityHelper.h>
#include <SelectionHelpers/SysReadSelectionHandle.h>
#include <SystematicsHandles/SysReadHandle.h>
#include <SystematicsHandles/SysWriteDecorHandle.h>
#include <SystematicsHandles/SysListHandle.h>
#include <SystematicsHandles/SysReadHandle.h>
#include <xAODJet/JetContainer.h>
#include <memory>

namespace CP
{
  /// \brief an algorithm for calling \ref IBTaggingEfficiencyTool

  class BTaggingEfficiencyAlg final : public EL::AnaAlgorithm
  {
    /// \brief the standard constructor
  public:
    BTaggingEfficiencyAlg (const std::string& name, 
                           ISvcLocator* pSvcLocator);


  public:
    StatusCode initialize () override;

  public:
    StatusCode execute () override;
    


    /// \brief the smearing tool
  private:
    ToolHandle<IBTaggingEfficiencyTool> m_efficiencyTool;

    /// \brief the systematics list we run
  private:
    SysListHandle m_systematicsList {this};

    /// \brief the jet collection we run on
  private:
    SysReadHandle<xAOD::JetContainer> m_jetHandle {
      this, "jets", "Jets", "the jet collection to run on"};

    /// \brief the preselection we apply to our input
  private:
    SysReadSelectionHandle m_preselection {
      this, "preselection", "", "the preselection to apply"};

    /// \brief the helper for OutOfValidity results
  private:
    OutOfValidityHelper m_outOfValidity {this};

    /// \brief the decoration for the b-tagging scale factor
  private:
    SysWriteDecorHandle<float> m_scaleFactorDecoration {
      this, "scaleFactorDecoration", "", "the decoration for the b-tagging efficiency scale factor"};

    /// \brief the decoration for the b-tagging selection
  private:
    SysReadSelectionHandle m_selectionHandle {
      this, "selectionDecoration", "", "the decoration for the asg selection"};

    /// \brief only run the inefficency for all jets
  private:
    bool m_onlyInefficiency {false};
  };
}

#endif
