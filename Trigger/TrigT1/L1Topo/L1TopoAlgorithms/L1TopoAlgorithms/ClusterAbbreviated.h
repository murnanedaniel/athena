//  ClusterAbbreviated.h
//  TopoCore
//  Created by Joerg Stelzer on 11/10/12.
//  Copyright (c) 2012 Joerg Stelzer. All rights reserved.

#ifndef TCS__ClusterAbbreviated
#define TCS__ClusterAbbreviated

#include "L1TopoInterfaces/SortingAlg.h"
#include "L1TopoEvent/TOBArray.h"

#include <iostream>
#include <vector>

namespace TCS {
   
   class ClusterAbbreviated : public SortingAlg {
   public:
      
      // constructor
      ClusterAbbreviated(const std::string & name);

      // destructor
      virtual ~ClusterAbbreviated();

      virtual TCS::StatusCode sort(const InputTOBArray & input, TOBArray & output);
      
   };

} // end of namespace TCS

#endif /* defined(__TopoCore__SortingAlg__) */
