/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/*! \file BlackBin.h file declares the dqm_algorithms::BlackBin  class.
 * \author Martin White
*/

#ifndef DQM_ALGORITHMS_BLACKBIN1D_H
#define DQM_ALGORITHMS_BLACKBIN1D_H

#include <dqm_core/Algorithm.h>

namespace dqm_algorithms
{
	struct BlackBin1D : public dqm_core::Algorithm
        {
	  BlackBin1D();

	    //overwrites virtual functions
	  BlackBin1D * clone( );
	  dqm_core::Result * execute( const std::string & , const TObject & , const dqm_core::AlgorithmConfig & );
          using dqm_core::Algorithm::printDescription;
	  void  printDescription(std::ostream& out);
 
	  private:
	  std::string name_;
	};
}

#endif // DQM_ALGORITHMS_BLACKBIN1D_H
