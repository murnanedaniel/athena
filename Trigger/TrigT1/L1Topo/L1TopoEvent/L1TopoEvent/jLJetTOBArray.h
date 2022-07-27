/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


#ifndef JLJETTOBARRAY_H
#define JLJETTOBARRAY_H

#include <iostream>
#include "L1TopoEvent/InputTOBArray.h"
#include "L1TopoEvent/DataArrayImpl.h"
#include "L1TopoEvent/jLJetTOB.h"
#include <vector>

namespace TCS {

   class TOBArray;

   class jLJetTOBArray : public InputTOBArray, public DataArrayImpl<jLJetTOB> {
   public:
      // default constructor
      jLJetTOBArray(const std::string & name, unsigned int reserve);

      virtual unsigned int size() const { return DataArrayImpl<jLJetTOB>::size(); }

      TOBArray asTOBArray() const;

   private:      
      // print method can be invoked via <<
      virtual void print(std::ostream&) const;
   };
   
}

#endif
