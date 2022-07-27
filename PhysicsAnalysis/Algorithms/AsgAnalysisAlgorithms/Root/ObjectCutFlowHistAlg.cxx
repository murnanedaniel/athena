/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack


//
// includes
//

#include <AsgAnalysisAlgorithms/ObjectCutFlowHistAlg.h>

#include <RootCoreUtils/StringUtil.h>
#include <TH1.h>

//
// method implementations
//

namespace CP
{
  ObjectCutFlowHistAlg ::
  ObjectCutFlowHistAlg (const std::string& name, 
                          ISvcLocator* pSvcLocator)
    : AnaAlgorithm (name, pSvcLocator)
  {
    declareProperty ("histPattern", m_histPattern, "the pattern for histogram names");

    declareProperty ("selection", m_selection, "the list of selection decorations");
    declareProperty ("selectionNCuts", m_selectionNCuts, "the number of cuts for each selection decoration");
  }



  StatusCode ObjectCutFlowHistAlg ::
  initialize ()
  {
    ANA_CHECK (m_inputHandle.initialize (m_systematicsList));
    ANA_CHECK (m_preselection.initialize (m_systematicsList, m_inputHandle, SG::AllowEmpty));
    ANA_CHECK (m_systematicsList.initialize());

    if (m_selectionNCuts.size() != m_selection.size())
    {
      ANA_MSG_ERROR ("selection and selectionNCuts properties need to be the same size");
      return StatusCode::FAILURE;
    }

    // Total label
    m_labels.push_back ("total");

    for (std::size_t iter = 0, end = m_selection.size(); iter != end; ++ iter)
    {
      const unsigned ncuts = m_selectionNCuts[iter];
      if (ncuts == 0)
      {
        ANA_MSG_ERROR ("all entries of selectionNCuts need to be greater than 0");
        return StatusCode::FAILURE;
      }
      if (ncuts > 8 * sizeof (SelectionType))
      {
        ANA_MSG_ERROR ("entries in selectionNCuts need to be less or equal to " << (8 * sizeof (SelectionType)));
        return StatusCode::FAILURE;
      }
      std::unique_ptr<ISelectionReadAccessor> accessor;
      ANA_CHECK (makeSelectionReadAccessor (m_selection[iter], accessor));
      m_accessors.push_back (std::make_pair (std::move (accessor), ncuts));
      for (unsigned i = 1; i <= ncuts; i++)
      {
        m_labels.push_back (m_accessors.back().first->label() + std::to_string(i));
      }
      m_allCutsNum += ncuts;
    }

    return StatusCode::SUCCESS;
  }



  StatusCode ObjectCutFlowHistAlg ::
  execute ()
  {
    for (const auto& sys : m_systematicsList.systematicsVector())
    {
      const xAOD::IParticleContainer *input = nullptr;
      ANA_CHECK (m_inputHandle.retrieve (input, sys));

      auto histIter = m_hist.find (sys);
      if (histIter == m_hist.end())
      {
        std::string name;
        ANA_CHECK (m_systematicsList.service().makeSystematicsName (name, m_histPattern, sys));

        ANA_CHECK (book (TH1F (name.c_str(), "object cut flow", m_allCutsNum+1, 0, m_allCutsNum+1)));

        m_hist.insert (std::make_pair (sys, hist (name)));
        histIter = m_hist.find (sys);
        assert (histIter != m_hist.end());

        for (unsigned i = 0; i < m_allCutsNum+1; i++)
        {
          histIter->second->GetXaxis()->SetBinLabel(i + 1, m_labels[i].c_str());
        }
      }

      for (const xAOD::IParticle *particle : *input)
      {
        if (m_preselection.getBool (*particle, sys))
        {
          bool keep = true;
          unsigned cutIndex = 1;
          histIter->second->Fill (0);
          for (const auto& accessor : m_accessors)
          {
            const auto selection = accessor.first->getBits (*particle);
            for (unsigned index = 0, end = accessor.second;
                 index != end; ++ index, ++ cutIndex)
            {
              if (selection & (1 << index))
              {
                histIter->second->Fill (cutIndex);
              } else
              {
                keep = false;
                break;
              }
            }
            if (!keep)
              break;
          }
        }
      }
    }

    return StatusCode::SUCCESS;
  }
}
