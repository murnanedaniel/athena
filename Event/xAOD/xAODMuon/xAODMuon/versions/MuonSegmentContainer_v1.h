/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef XAODMUON_VERSIONS_MUONSEGMENTCONTAINER_V1_H
#define XAODMUON_VERSIONS_MUONSEGMENTCONTAINER_V1_H

// Core include(s):
#include "AthContainers/DataVector.h"

// Local include(s):
#include "xAODMuon/versions/MuonSegment_v1.h"

namespace xAOD {
   /// The container is a simple typedef for now
   typedef DataVector< xAOD::MuonSegment_v1 > MuonSegmentContainer_v1;
}

// Set up a CLID for the container:
#include "xAODCore/CLASS_DEF.h"
CLASS_DEF( xAOD::MuonSegmentContainer_v1, 1129401482, 1 )

#endif // XAODMUON_VERSIONS_MUONSEGMENTCONTAINER_V1_H
