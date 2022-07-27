/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DQM_ALGORITHMS_CORRELATIONYX_H
#define DQM_ALGORITHMS_CORRELATIONYX_H



#include "dqm_core/Algorithm.h"
#include <string>
#include <iosfwd>

namespace dqm_algorithms {

  class CorrelationYX : public dqm_core::Algorithm {
  public:

    CorrelationYX();
  
    virtual ~CorrelationYX();
    virtual dqm_core::Algorithm*  clone();
    virtual dqm_core::Result*     execute( const std::string& name, const TObject& object,
					   const dqm_core::AlgorithmConfig& config );
    using dqm_core::Algorithm::printDescription;
    virtual void                  printDescription(std::ostream& out);

  private:
    std::string  m_name;
  };

} //namespace dqm_algorithms

#endif // DQM_ALGORITHMS_CORRELATIONYX_H
