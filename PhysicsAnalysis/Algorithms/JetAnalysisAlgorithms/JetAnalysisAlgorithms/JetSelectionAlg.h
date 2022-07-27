/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack


#ifndef JET_ANALYSIS_ALGORITHMS__JET_SELECTION_ALG_H
#define JET_ANALYSIS_ALGORITHMS__JET_SELECTION_ALG_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <JetInterface/IJetSelector.h>
#include <SelectionHelpers/SysReadSelectionHandle.h>
#include <SelectionHelpers/SysWriteSelectionHandle.h>
#include <SystematicsHandles/SysListHandle.h>
#include <SystematicsHandles/SysReadHandle.h>
#include <xAODJet/JetContainer.h>

namespace CP
{
  /// \brief an algorithm for calling \ref IJetSelector

  class JetSelectionAlg final : public EL::AnaAlgorithm
  {
    /// \brief the standard constructor
  public:
    JetSelectionAlg (const std::string& name, 
                         ISvcLocator* pSvcLocator);


  public:
    StatusCode initialize () override;

  public:
    StatusCode execute () override;
    


    /// \brief the selection tool
  private:
    ToolHandle<IJetSelector> m_selectionTool;

    /// \brief the systematics list we run
  private:
    SysListHandle m_systematicsList {this};

    /// \brief the jet collection we run on
  private:
    SysReadHandle<xAOD::JetContainer> m_jetHandle {
      this, "jets", "AntiKt4EMTopoJets", "the jet collection to run on"};

    /// \brief the preselection we apply to our input
  private:
    SysReadSelectionHandle m_preselection {
      this, "preselection", "", "the preselection to apply"};

    /// \brief the decoration for the jet selection
  private:
    SysWriteSelectionHandle m_selectionHandle {
      this, "selectionDecoration", "clean_jet", "the decoration for the jet selection"};
  };
}

#endif
