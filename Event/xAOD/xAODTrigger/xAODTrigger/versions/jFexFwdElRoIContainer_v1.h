/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef XAODTRIGGER_VERSIONS_JFEXFWDELROICONTAINER_V1_H
#define XAODTRIGGER_VERSIONS_JFEXFWDELROICONTAINER_V1_H

// Core include(s):
#include "AthContainers/DataVector.h"

// Local include(s):
#include "xAODTrigger/versions/jFexFwdElRoI_v1.h"

 namespace xAOD {
    // Declare the Jet TOB container type
       typedef DataVector< xAOD::jFexFwdElRoI_v1 > jFexFwdElRoIContainer_v1;
       }

#endif // XAODTRIGGER_VERSIONS_JFEXFWDELROICONTAINER_V1_H
