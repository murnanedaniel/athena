///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MCPARTICLEEVENTTPCNV_ROOTTRUTHPARTICLECNVTOOL_H
#define MCPARTICLEEVENTTPCNV_ROOTTRUTHPARTICLECNVTOOL_H

#include "McParticleKernel/ITruthParticleCnvTool.h"

class RootTruthParticleCnvTool
  : public ITruthParticleCnvTool
{
public:
  /** Converts a @c McEventCollection into an @c TruthParticleContainer (ie:
   *  converts it into an AOD compliant collection).
   *  @in  mcEvts the @c McEventCollection holding the @c HepMC::GenEvent we
   *      want to convert into a @c TruthParticleContainer
   *  @in  genEvtIndex the index to the @c HepMC::GenEvent to be converted
   *  @out mcParts a valid pointer to a @c TruthParticleContainer which will
   *       be filled with adaptors to @c HepMC::GenParticles.
   */
  StatusCode convert( const McEventCollection* mcEvts,
		      const unsigned int genEvtIndex,
		      TruthParticleContainer* mcParts,
		      const ITruthParticleVisitor* visitor ) const;

  /** Helper method to get the charge of a particle given its PDG Id.
   */
  double chargeFromPdgId( int pdgId ) const;

  virtual StatusCode queryInterface(const InterfaceID& riid,
                                    void** ppvInterface);
  virtual unsigned long addRef();
  virtual unsigned long release();
  virtual StatusCode setProperty( const Property& p );
  virtual StatusCode setProperty( const std::string& s );
  virtual StatusCode setProperty( const std::string& n, const std::string& v );
  virtual StatusCode getProperty( Property* p ) const;
  virtual const Property& getProperty( const std::string& name) const;
  virtual StatusCode getProperty( const std::string& n, std::string& v ) const;
  virtual const std::vector<Property*>& getProperties( ) const;
  virtual bool hasProperty(const std::string& name) const;

  virtual const std::string&  type() const;
  virtual const IInterface*   parent() const;
  virtual StatusCode configure();
  virtual StatusCode initialize();
  virtual StatusCode sysInitialize();
  virtual StatusCode reinitialize();
  virtual StatusCode sysReinitialize();
  virtual StatusCode start();
  virtual StatusCode sysStart();
  virtual StatusCode restart();
  virtual StatusCode sysRestart();
  virtual StatusCode stop();
  virtual StatusCode sysStop();
  virtual StatusCode finalize();
  virtual StatusCode sysFinalize();
  virtual StatusCode terminate();
  virtual unsigned long refCount() const;
  virtual const std::string& name() const;
  virtual StatusCode execute();
#ifdef GAUDIKERNEL_STATEMACHINE_H_
  virtual Gaudi::StateMachine::State FSMState() const;
#endif

#ifdef ATHENAHIVE
  virtual const DataObjectDescriptorCollection & inputDataObjects() const;
  virtual const DataObjectDescriptorCollection & outputDataObjects() const;
#endif


};

#endif // not MCPARTICLEEVENTTPCNV_ROOTTRUTHPARTICLECNVTOOL_H
