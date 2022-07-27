/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Tadej Novak

#ifndef ASG_ANALYSIS_ALGORITHMS__EVENT_FLAG_SELECTION_ALG_H
#define ASG_ANALYSIS_ALGORITHMS__EVENT_FLAG_SELECTION_ALG_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <EventBookkeeperTools/FilterReporterParams.h>
#include <SelectionHelpers/ISelectionReadAccessor.h>

namespace CP
{
  /// \brief an algorithm for selecting events by flag
  class EventFlagSelectionAlg final : public EL::AnaAlgorithm
  {
  public:
    EventFlagSelectionAlg(const std::string &name,
                          ISvcLocator *svcLoc = nullptr);

    virtual StatusCode initialize() final;
    virtual StatusCode execute() final;
    virtual StatusCode finalize() final;

  private:
    /// \brief flags that we want to select events with
    std::vector<std::string> m_selFlags;
    
    /// \brief invert flags
    std::vector<bool> m_invertFlags;
    
    /// \brief a vector of accessors to read the flags
    std::vector<std::unique_ptr<ISelectionReadAccessor>> m_accessors;

    /// \brief the filter reporter params
    FilterReporterParams m_filterParams {this, "EventFlagSelection", "event flag selection"};
  };
}

#endif
