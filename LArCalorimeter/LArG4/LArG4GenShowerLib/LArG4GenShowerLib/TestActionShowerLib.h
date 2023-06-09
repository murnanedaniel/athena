/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


#ifndef TestActionShowerLib_H
#define TestActionShowerLib_H

// STL includes
#include <string>
#include <vector>

// athena simulation includes

#include "G4AtlasTools/UserActionBase.h"

// forward declarations in namespaces
namespace ShowerLib {
  class StepInfoCollection;
}
namespace HepMC {
  class GenParticle;
}
// forward declarations in global namespace
//class StoreGateSvc;
class LArVCalculator;
class G4VSolid;
class G4AffineTransform;

  /**
   *
   *   @short Class for collecting G4 hit information
   *
   *          Collect and store Geant4 hit information, i.e.
   *          position, deposited energy and time, from hits
   *          for the creation of a shower library
   *
   *  @author Wolfgang Ehrenfeld, University of Hamburg, Germany
   *  @author Sasha Glazov, DESY Hamburg, Germany
   *
   * @version \$Id: TestActionShowerLib.h 746605 2016-05-12 13:25:40Z disimone $
   *
   */

class TestActionShowerLib final: public UserActionBase {

 public:

  //! default constructor
  TestActionShowerLib(const std::string& type, const std::string& name, const IInterface* parent);

  ~TestActionShowerLib();

  //! run code at begin of event
  void BeginOfEvent(const G4Event*) override;
  //! run code at end of event
  void EndOfEvent(const G4Event*) override;
  //! run code at begin of run
  void BeginOfRun(const G4Run*) override;
  //! run code at end of event
  void EndOfRun(const G4Run*) override;
  //! run code after each step
  void Step(const G4Step*) override;

  virtual StatusCode queryInterface(const InterfaceID&, void**) override;

 private:

  /* data members */

  LArVCalculator* m_current_calculator;
  G4VSolid* m_current_solid;
  G4AffineTransform* m_current_transform;

  // calculators 
  LArVCalculator* m_calculator_EMECIW;            //!< pointer to EMEC inner wheel calculator
  LArVCalculator* m_calculator_EMECOW;            //!< pointer to EMEC outer wheel calculator
  

  ShowerLib::StepInfoCollection* m_eventSteps;    //!< collection of StepInfo

};


#include "G4AtlasInterfaces/IBeginEventAction.h"
#include "G4AtlasInterfaces/IEndEventAction.h"
#include "G4AtlasInterfaces/IBeginRunAction.h"
#include "G4AtlasInterfaces/IEndRunAction.h"
#include "G4AtlasInterfaces/ISteppingAction.h"

#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/ServiceHandle.h"
namespace G4UA{
  
  
  class TestActionShowerLib:
  public IBeginEventAction,  public IEndEventAction,  public IBeginRunAction,  public IEndRunAction,  public ISteppingAction
  {
    
  public:
    TestActionShowerLib();
    virtual void beginOfEvent(const G4Event*) override;
    virtual void endOfEvent(const G4Event*) override;
    virtual void beginOfRun(const G4Run*) override;
    virtual void endOfRun(const G4Run*) override;
    virtual void processStep(const G4Step*) override;
  private:
    
    typedef ServiceHandle<StoreGateSvc> StoreGateSvc_t;
    /// Pointer to StoreGate (event store by default)
    mutable StoreGateSvc_t m_evtStore;
    /// Pointer to StoreGate (detector store by default)
    mutable StoreGateSvc_t m_detStore;
    
    /* data members */
    
    LArVCalculator* m_current_calculator;
    G4VSolid* m_current_solid;
    G4AffineTransform* m_current_transform;
    
    // calculators 
    LArVCalculator* m_calculator_EMECIW;            //!< pointer to EMEC inner wheel calculator
    LArVCalculator* m_calculator_EMECOW;            //!< pointer to EMEC outer wheel calculator
    
    
    ShowerLib::StepInfoCollection* m_eventSteps;    //!< collection of StepInfo

}; // class TestActionShowerLib


} // namespace G4UA 


#endif // TestActionShowerLib_H
