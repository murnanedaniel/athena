/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// IJetModifier.h

// David Adams
// January 2014

/// IJetModifier is a dual-use tool interface for a tool that
/// modifies a jet collection.

#ifndef IJetModifier_H
#define IJetModifier_H

#include "AsgTools/IAsgTool.h"
#include "xAODJet/JetContainer.h"

class IJetModifier : virtual public asg::IAsgTool {
ASG_TOOL_INTERFACE(IJetModifier)

public:

  /// Destructor.
  virtual ~IJetModifier() { };

  /// Method to modify a jet collection.
  /// Returns 0 for success.
  virtual int modify(xAOD::JetContainer& jets) const =0;

  /// Method to return the list of input containers.
  /// The names of required input containers are appended to connames.
  /// Returns nonzero for error.
  /// Default returns 0 and adds no names.
  virtual int inputContainerNames(std::vector<std::string>& connames);

};

#endif
