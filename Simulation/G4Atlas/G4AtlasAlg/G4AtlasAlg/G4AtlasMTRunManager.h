/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef G4ATLASALG_G4ATLASMTRUNMANAGER_H
#define G4ATLASALG_G4ATLASMTRUNMANAGER_H

// Hide multi-threading classes from builds without G4MT
#include "G4Types.hh"
#ifdef G4MULTITHREADED

// Geant4 includes
#include "G4MTRunManager.hh"

// Framework includes
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ServiceHandle.h"
#include "AthenaBaseComps/AthMessaging.h"

// G4Atlas includes
#include "G4AtlasInterfaces/IPhysicsListSvc.h"
#include "G4AtlasInterfaces/IDetectorGeometrySvc.h"
#include "G4AtlasInterfaces/IFastSimulationMasterTool.h"


/// @class G4AtlasMTRunManager
/// @brief ATLAS master run manager for master-slave multi-threading model
///
/// This thread-local singleton (eww) is used on the G4MT master thread
/// to do some setup for the G4 run.
///
/// The corresponding worker thread run manager is G4AtlasWorkerRunManager.
///
/// @author Steve Farrell <Steven.Farrell@cern.ch>
///
class G4AtlasMTRunManager: public G4MTRunManager, public AthMessaging {

public:

  /// Get the (pure) singleton instance
  static G4AtlasMTRunManager* GetG4AtlasMTRunManager();

  /// G4 function called at end of run
  void RunTermination() override final;

  /// We cram all of the initialization of the run manager stuff in here.
  /// This then includes some of the things that in normal G4 are called
  /// immediately before the event loop.
  void Initialize() override final;

  /// Disable G4's barrier synchronization by implementing these methods
  /// and leaving them empty
  virtual void ThisWorkerReady() override final {};
  virtual void ThisWorkerEndEventLoop() override final {};

  /// Configure the detector geometry service handle
  void SetDetGeoSvc(const std::string& typeAndName) {
    m_detGeoSvc.setTypeAndName(typeAndName);
  }

  /// Configure the Fast Simulation Master Tool handle
  void SetFastSimMasterTool(const std::string& typeAndName) {
    m_fastSimTool.setTypeAndName(typeAndName);
  }

  /// Configure the Physics List Tool handle
  void SetPhysListSvc(const std::string& typeAndName) {
    m_physListSvc.setTypeAndName(typeAndName);
  }

 protected:

  /// Initialize the G4 geometry on the master
  void InitializeGeometry() override final;

  // Initialize the physics list on the master
  void InitializePhysics() override final;

  /// Disable G4's barrier synchronization by implementing these methods
  /// and leaving them empty
  virtual void WaitForReadyWorkers() override final {};
  virtual void WaitForEndEventLoopWorkers() override final {};

private:

  /// Pure singleton private constructor
  G4AtlasMTRunManager();

private:
  /// Handle to the detector geometry service.
  /// Not ideal, because we can't configure this.
  ServiceHandle<IDetectorGeometrySvc> m_detGeoSvc;

  /// Handle to the physics list tool.
  /// Not ideal, because we can't configure this.
  ServiceHandle<IPhysicsListSvc> m_physListSvc;

  /// Handle to the fast sim tool.
  /// Not ideal, because we can't configure this.
  ToolHandle<IFastSimulationMasterTool> m_fastSimTool;

}; // class G4AtlasMTRunManager

#endif // G4MULTITHREADED

#endif
