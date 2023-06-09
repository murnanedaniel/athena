//  GenericListGenerator.h
//  TopoCore
//  Created by Joerg Stelzer on 11/10/12.
//  Copyright (c) 2012 Joerg Stelzer. All rights reserved.

#ifndef TCS__GenericListGenerator
#define TCS__GenericListGenerator

#include "L1TopoInterfaces/SortingAlg.h"
#include "L1TopoEvent/TOBArray.h"

#include <iostream>
#include <vector>

namespace TCS {
   
   class GenericListGenerator : public SortingAlg {
   public:
      
      // constructor
      GenericListGenerator(const std::string & name);

      // destructor
      virtual ~GenericListGenerator();

      virtual TCS::StatusCode sort(const InputTOBArray & input, TOBArray & output);

   protected:
    
      JetTOB::JetSize m_jetsize;
      
   };

} // end of namespace TCS

#endif /* defined(__TopoCore__SortingAlg__) */
