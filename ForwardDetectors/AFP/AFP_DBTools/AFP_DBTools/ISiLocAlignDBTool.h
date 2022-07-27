/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/// @brief  Header file for interface of SiLocAlignDBTool used to read local alignment for database
#ifndef AFP_DBTOOLS_ISILOCALIGNDBTOOL_H
#define AFP_DBTOOLS_ISILOCALIGNDBTOOL_H


// FrameWork includes
#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/EventContext.h"

// forward declaration of nlohmann::json
#include "nlohmann/json_fwd.hpp"

namespace AFP
{
  // forward declarations
  class SiLocAlignData;
  
  /// Interface to tool providing local alignment of silicon detectors from the conditions database.
  class ISiLocAlignDBTool : virtual public IAlgTool
  {
  public:
    DeclareInterfaceID(ISiLocAlignDBTool, 1, 0);

    /// Provide alignment parameters for a given plane. Returns zeros if no data available.    
    virtual nlohmann::json alignmentData(const EventContext& ctx) const = 0;
    virtual const SiLocAlignData alignment(const nlohmann::json& jsondata, const int stationID, const int planeID) const = 0;
  };

}      // namespace AFP

#endif // > ! AFP_DBTOOLS_ISILOCALIGNDBTOOL_H
