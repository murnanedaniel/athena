/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "PhotonKillerTool.h"

namespace G4UA
{

  //---------------------------------------------------------------------------
  PhotonKillerTool::PhotonKillerTool(const std::string& type,
                                     const std::string& name,
                                     const IInterface* parent)
    : UserActionToolBase<PhotonKiller>(type, name, parent)
  {
  }

  //---------------------------------------------------------------------------
  std::unique_ptr<PhotonKiller>
  PhotonKillerTool::makeAndFillAction(G4AtlasUserActions& actionList)
  {
    ATH_MSG_DEBUG("Making a PhotonKiller action");
    auto action = std::make_unique<PhotonKiller>();
    actionList.trackingActions.push_back( action.get() );
    actionList.steppingActions.push_back( action.get() );
    return action;
  }

} // namespace G4UA
