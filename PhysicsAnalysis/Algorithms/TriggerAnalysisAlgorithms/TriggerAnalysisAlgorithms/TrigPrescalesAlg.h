/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Tadej Novak


#ifndef TRIGGER_ANALYSIS_ALGORITHMS__TRIG_PRESCALES_ALG_H
#define TRIGGER_ANALYSIS_ALGORITHMS__TRIG_PRESCALES_ALG_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <AsgAnalysisInterfaces/IPileupReweightingTool.h>

namespace CP
{
  /// \brief an algorithm for retrieving trigger prescales

  class TrigPrescalesAlg final : public EL::AnaAlgorithm
  {
    /// \brief the standard constructor
  public:
    TrigPrescalesAlg (const std::string& name, 
                      ISvcLocator* pSvcLocator);


  public:
    StatusCode initialize () override;

  public:
    StatusCode execute () override;


    /// \brief the pile-up reweighting tool
  private:
    ToolHandle<IPileupReweightingTool> m_pileupReweightingTool;

    /// \brief list of prescaled triggers or trigger chains
  private:
    std::vector<std::string> m_trigList;
  
    /// \brief list of all triggers or trigger chains
  private:
    std::vector<std::string> m_trigListAll;

    /// \brief the decoration for trigger prescales
  private:
    std::string m_prescaleDecoration;

    /// \brief the accessors for \ref m_prescaleDecoration and \ref m_trigList combination
  private:
    std::vector<SG::AuxElement::Decorator<float>> m_prescaleAccessors;
  };
}

#endif
