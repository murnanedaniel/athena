/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef _MYDATAOBJ_
#define _MYDATAOBJ_

//sample data class for the Read/Write example
//it just wraps an int. Notice that is does not inherit from Gaudi DataObject

//a dummy base class to test the symlinks
class BaseClass { };

class MyDataObj : public BaseClass {
 
public: 
 
  MyDataObj(): BaseClass(), m_val(0) {};
  MyDataObj(int i): m_val(i) {};
  virtual ~MyDataObj(){};

  void val(int i) { m_val = i; }
  int val() const { return m_val; }

private:
  int m_val;
};
//using the macros below we can assign an identifier (and a version) 
//to the type MyDataObj
//This is required and checked at compile time when you try to record/retrieve
#include "SGTools/CLASS_DEF.h"
CLASS_DEF(MyDataObj, 8000, 1)
CLASS_DEF(BaseClass, 1434, 1)

#endif


