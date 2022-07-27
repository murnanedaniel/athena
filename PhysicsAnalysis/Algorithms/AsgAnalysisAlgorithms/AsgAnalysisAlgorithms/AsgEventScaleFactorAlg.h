/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/// @author Tadej Novak


#ifndef ASG_ANALYSIS_ALGORITHMS__ASG_EVENT_SCALE_FACTOR_ALG_H
#define ASG_ANALYSIS_ALGORITHMS__ASG_EVENT_SCALE_FACTOR_ALG_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <SelectionHelpers/SysReadSelectionHandle.h>
#include <SystematicsHandles/SysReadHandle.h>
#include <SystematicsHandles/SysReadDecorHandle.h>
#include <SystematicsHandles/SysWriteDecorHandle.h>
#include <SystematicsHandles/SysListHandle.h>
#include <xAODBase/IParticleContainer.h>
#include <xAODEventInfo/EventInfo.h>

namespace CP
{
  /// \brief an algorithm for calculating per-event scale factors

  class AsgEventScaleFactorAlg final : public EL::AnaAlgorithm
  {
    /// \brief the standard constructor
  public:
    AsgEventScaleFactorAlg (const std::string& name, 
                            ISvcLocator* pSvcLocator);


  public:
    StatusCode initialize () override;

  public:
    StatusCode execute () override;


    /// \brief the systematics list we run
  private:
    SysListHandle m_systematicsList {this};

    /// \brief the event info we run on (empty by default)
  private:
    SysReadHandle<xAOD::EventInfo> m_eventInfoHandle {
      this, "eventInfo", "EventInfo", "the event info object to run on"};

    /// \brief the jet collection we run on
  private:
    SysReadHandle<xAOD::IParticleContainer> m_particleHandle {
      this, "particles", "", "the particle collection to run on"};

    /// \brief the preselection we apply to our input
  private:
    SysReadSelectionHandle m_preselection {
      this, "preselection", "", "the preselection to apply"};

    /// \brief the decoration for reading the scale factor
  private:
    SysReadDecorHandle<float> m_scaleFactorInputDecoration {
      this, "scaleFactorInputDecoration", "", "the decoration for the input efficiency scale factor"};

    /// \brief the decoration for writing the scale factor
  private:
    SysWriteDecorHandle<float> m_scaleFactorOutputDecoration {
      this, "scaleFactorOutputDecoration", "", "the decoration for the output efficiency scale factor"};
  };
}

#endif
