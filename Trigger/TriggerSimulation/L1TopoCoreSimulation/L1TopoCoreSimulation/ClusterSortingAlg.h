//
//  ClusterSortingAlg.h
//  TopoCore
//
//  Created by Joerg Stelzer on 11/10/12.
//  Copyright (c) 2012 Joerg Stelzer. All rights reserved.
//

#ifndef __TopoCore__ClusterSortingAlg__
#define __TopoCore__ClusterSortingAlg__

#include "L1TopoCoreSimulation/SortingAlg.h"
#include "L1TopoCoreSimulation/TOBArray.h"

#include <iostream>
#include <vector>

namespace TCS {
   
   class ClusterSortingAlg : public SortingAlg {
   public:
      
      // constructor
      ClusterSortingAlg(const std::string & name);

      // destructor
      virtual ~ClusterSortingAlg(){};

      virtual TCS::StatusCode sort(const InputTOBArray & input, TOBArray & output);
      
   };

} // end of namespace TCS

#endif /* defined(__TopoCore__SortingAlg__) */
