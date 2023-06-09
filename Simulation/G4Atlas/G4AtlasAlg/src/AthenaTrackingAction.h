/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef G4AtlasAlg_AthenaTrackingAction_H
#define G4AtlasAlg_AthenaTrackingAction_H

#include "G4AtlasTools/UserActionBase.h"

/// @class AthenaTrackingAction
/// @brief User action for pre/post tracking truth handling.
///
/// This is the current implementation based on the (V1) user action design
/// as part of the simulation infrastructure migrations. The multi-threaded
/// (V2) design is still in the works.
///
class AthenaTrackingAction : public UserActionBase {
public:

  /// Standard constructor
  AthenaTrackingAction(const std::string& type, const std::string& name,
                       const IInterface* parent);

  /// Geant4 method called at the beginning of tracking a particle.
  virtual void PreTracking(const G4Track*) override;
  /// Geant4 method called at the end of tracking a particle.
  virtual void PostTracking(const G4Track*) override;

  /// Gaudi boiler-plate.
  virtual StatusCode queryInterface(const InterfaceID&, void**) override;

}; // class AthenaTrackingAction


//=============================================================================
// New design for multi-threading (V2 migration) follows.
//=============================================================================

#include "AthenaKernel/MsgStreamMember.h"
#include "G4AtlasInterfaces/IPreTrackingAction.h"
#include "G4AtlasInterfaces/IPostTrackingAction.h"

namespace G4UA
{

  /// @class AthenaTrackingAction
  /// @brief User action for pre/post tracking truth handling.
  ///
  /// This is the new/upcoming implementation for multi-threaded simulation.
  ///
  class AthenaTrackingAction : public IPreTrackingAction,
                               public IPostTrackingAction
  {

    public:

      /// Constructor
      AthenaTrackingAction(MSG::Level lvl = MSG::INFO);

      /// @brief Called before tracking a new particle.
      ///
      /// If the track meets certain conditions, we save it in the
      /// EventInformation and possibly construct a new AtlasTrajectory
      /// which will be used for writing out truth particles later.
      virtual void preTracking(const G4Track*) override;

      /// @brief Called after tracking a particle.
      ///
      /// Here we reset the AtlasTrajectory if it was created.
      virtual void postTracking(const G4Track*) override;

    private:

      /// Log a message using the Athena controlled logging system
      MsgStream& msg( MSG::Level lvl ) const { return m_msg << lvl; }
      bool msgLvl( MSG::Level lvl ) const { return m_msg.get().level() <= lvl; }
      mutable Athena::MsgStreamMember m_msg;

  }; // class AthenaTrackingAction

} // namespace G4UA

#endif
