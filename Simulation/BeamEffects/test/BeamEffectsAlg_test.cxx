/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @author John Chapman
 * @date June, 2016
 * @brief Tests for BeamEffectsAlg.
 */

#undef NDEBUG

// Framework
#include "TestTools/initGaudi.h"

// Google Test
#include "gtest/gtest.h"

// HepMC includes
#include "AtlasHepMC/GenEvent.h"
#include "AtlasHepMC/GenVertex.h"
#include "AtlasHepMC/GenParticle.h"
#include "AtlasHepMC/Operators.h"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

// Athena headers
#include "CxxUtils/checker_macros.h"
#include "GeneratorObjects/McEventCollection.h"

// Tested AthAlgorithm
#include "../src/BeamEffectsAlg.h"


namespace SimTesting {

  // needed every time an AthAlgorithm, AthAlgTool or AthService is instantiated
  ISvcLocator* g_svcLoc ATLAS_THREAD_SAFE = nullptr;

  // global test environment takes care of setting up Gaudi
  class GaudiEnvironment : public ::testing::Environment {
  protected:
    virtual void SetUp() override {
      Athena_test::initGaudi(SimTesting::g_svcLoc);
    }
  };

  class BeamEffectsAlg_test : public ::testing::Test {

  protected:
    virtual void SetUp() override {
      m_alg = new Simulation::BeamEffectsAlg{"BeamEffectsAlg", g_svcLoc};
      ASSERT_TRUE( g_svcLoc->service("StoreGateSvc", m_sg) );
    }

    virtual void TearDown() override {
      ASSERT_TRUE( m_alg->finalize().isSuccess() );
      delete m_alg;
    }

    //
    // accessors for private methods
    // NB: This works because BeamEffectsAlg_test is a friend of the tested
    //     BeamEffectsAlg AthAlgorithm
    //
    template<typename... Args>
    StatusCode patchSignalProcessVertex(Args&&... args) const {
      return m_alg->patchSignalProcessVertex(std::forward<Args>(args)...);
    }

    Simulation::BeamEffectsAlg* m_alg{};
    StoreGateSvc* m_sg{};
  };   // BeamEffectsAlg_test fixture


  TEST_F(BeamEffectsAlg_test, empty_alg_execute) {
    ASSERT_TRUE( m_alg->initialize().isSuccess() );
    EventContext ctx(0,0);
    ctx.setExtension( Atlas::ExtendedEventContext( m_sg, 0 ) );
    // expected to fail as input collection doesn't exist
    ASSERT_TRUE( m_alg->execute(ctx).isFailure() );
  }

  TEST_F(BeamEffectsAlg_test, set_properties) {
    // ordering A, C, B is on purpose to test for unintended alphabetic ordering
    std::string  inputPropertyValue = "'TestGEN_EVENT'";
    std::string outputPropertyValue = "'TestBeamTruthEvent'";
    ASSERT_TRUE( m_alg->setProperty( "InputMcEventCollection",   inputPropertyValue).isSuccess() );
    ASSERT_TRUE( m_alg->setProperty( "OutputMcEventCollection", outputPropertyValue).isSuccess() );
    ASSERT_TRUE( m_alg->setProperty( "ISFRun", true).isSuccess()  );
  }

  TEST_F(BeamEffectsAlg_test, patchSignalProcessVertex_empty_GenEvent) {
    HepMC::GenEvent ge;
    ASSERT_TRUE( patchSignalProcessVertex(ge).isSuccess() );
    ASSERT_TRUE( HepMC::signal_process_vertex(&ge)==nullptr );
  }

  TEST_F(BeamEffectsAlg_test, signal_process_vertex_exists) {
    HepMC::GenEvent ge;
    CLHEP::HepLorentzVector myPos( 1.0, 1.0, 1.0, 1.0);
    HepMC::GenVertexPtr  myVertex = HepMC::newGenVertexPtr( HepMC::FourVector(myPos.x(),myPos.y(),myPos.z(),myPos.t()), -1 );
    HepMC::set_signal_process_vertex(&ge, myVertex );
    ASSERT_TRUE( patchSignalProcessVertex(ge).isSuccess() );
    ASSERT_TRUE( HepMC::signal_process_vertex(&ge)==myVertex );
  }

  TEST_F(BeamEffectsAlg_test, add_signal_process_vertex_atlasG4) {
    HepMC::GenEvent ge;
    CLHEP::HepLorentzVector myPos( 1.0, 1.0, 1.0, 1.0);
    HepMC::GenVertexPtr  myVertex = HepMC::newGenVertexPtr( HepMC::FourVector(myPos.x(),myPos.y(),myPos.z(),myPos.t()), -1 );
    ge.add_vertex( myVertex );
    ASSERT_TRUE( HepMC::signal_process_vertex(&ge)==nullptr );
    ASSERT_TRUE( m_alg->setProperty( "ISFRun", false).isSuccess()  );
    ASSERT_TRUE( patchSignalProcessVertex(ge).isSuccess() );
    ASSERT_TRUE( HepMC::signal_process_vertex(&ge)==myVertex );
#ifdef HEPMC3
//Not needed for HepMC3
#else
    ASSERT_EQ( *HepMC::signal_process_vertex(&ge), *myVertex );
#endif
  }

  TEST_F(BeamEffectsAlg_test, add_signal_process_vertex_isfG4) {
    HepMC::GenEvent ge;
    CLHEP::HepLorentzVector myPos( 1.0, 1.0, 1.0, 1.0);
    HepMC::GenVertexPtr  myVertex = HepMC::newGenVertexPtr( HepMC::FourVector(myPos.x(),myPos.y(),myPos.z(),myPos.t()), -1 );
    HepMC::GenVertexPtr  dummyVertex = HepMC::newGenVertexPtr();
    ge.add_vertex( myVertex );
    ASSERT_TRUE( HepMC::signal_process_vertex(&ge)==nullptr );
    ASSERT_TRUE( m_alg->setProperty( "ISFRun", true).isSuccess()  );
    ASSERT_TRUE( patchSignalProcessVertex(ge).isSuccess() );
    ASSERT_TRUE( HepMC::signal_process_vertex(&ge)!=myVertex );
#ifdef HEPMC3
//Not needed for HepMC3
#else
    ASSERT_EQ( *HepMC::signal_process_vertex(&ge), *dummyVertex );
#endif
  }

  TEST_F(BeamEffectsAlg_test, execute_pass_through) {
    // create dummy input McEventCollection containing a dummy GenEvent
    EventContext ctx(0,0);
    ctx.setExtension( Atlas::ExtendedEventContext( m_sg, 0 ) );
    SG::WriteHandleKey<McEventCollection> inputTestDataKey{"GEN_EVENT"};
    ASSERT_TRUE( inputTestDataKey.initialize().isSuccess() );
    SG::WriteHandle<McEventCollection> inputTestDataHandle{inputTestDataKey, ctx};
    inputTestDataHandle = std::make_unique<McEventCollection>();
    inputTestDataHandle->push_back(new HepMC::GenEvent());
    HepMC::GenEvent& ge = *(inputTestDataHandle->at(0));
    CLHEP::HepLorentzVector myPos( 0.0, 0.0, 0.0, 0.0);
    HepMC::GenVertexPtr  myVertex = HepMC::newGenVertexPtr( HepMC::FourVector(myPos.x(),myPos.y(),myPos.z(),myPos.t()), -1 );
    HepMC::FourVector fourMomentum1( 0.0, 0.0, 1.0, 1.0*CLHEP::TeV);
    HepMC::GenParticlePtr  inParticle1 = HepMC::newGenParticlePtr(fourMomentum1, 2, 10);
    myVertex->add_particle_in(inParticle1);
    HepMC::FourVector fourMomentum2( 0.0, 0.0, -1.0, 1.0*CLHEP::TeV);
    HepMC::GenParticlePtr  inParticle2 = HepMC::newGenParticlePtr(fourMomentum2, -2, 10);
    myVertex->add_particle_in(inParticle2);
    HepMC::FourVector fourMomentum3( 0.0, 1.0, 0.0, 1.0*CLHEP::TeV);
    HepMC::GenParticlePtr  inParticle3 = HepMC::newGenParticlePtr(fourMomentum3, 2, 10);
    myVertex->add_particle_out(inParticle3);
    HepMC::FourVector fourMomentum4( 0.0, -1.0, 0.0, 1.0*CLHEP::TeV);
    HepMC::GenParticlePtr  inParticle4 = HepMC::newGenParticlePtr(fourMomentum4, -2, 10);
    myVertex->add_particle_out(inParticle4);
    ge.add_vertex( myVertex );
    HepMC::set_signal_process_vertex(&ge, myVertex );
    ge.set_beam_particles(inParticle1,inParticle2);
    //
    ASSERT_TRUE( m_alg->initialize().isSuccess() );
    ASSERT_TRUE( m_alg->execute(ctx).isSuccess() );
    SG::ReadHandleKey<McEventCollection>     outputTestDataKey{"BeamTruthEvent"};
    ASSERT_TRUE( outputTestDataKey.initialize().isSuccess() );
    SG::ReadHandle<McEventCollection>     outputTestDataHandle{outputTestDataKey,ctx};
    ASSERT_TRUE( outputTestDataHandle.isValid() );
#ifdef HEPMC3
//This should compare the content 
    ASSERT_EQ(*(HepMC::signal_process_vertex(outputTestDataHandle->at(0))), *(HepMC::signal_process_vertex((const HepMC::GenEvent*)inputTestDataHandle->at(0))));
    ASSERT_EQ(**(outputTestDataHandle->at(0)->vertices().begin()), **(((const HepMC::GenEvent*)inputTestDataHandle->at(0))->vertices().begin()));
    ASSERT_EQ(*(outputTestDataHandle->at(0)->beams().at(0)),*(inputTestDataHandle->at(0)->beams().at(0)));
    ASSERT_EQ(*(outputTestDataHandle->at(0)->beams().at(1)),*(inputTestDataHandle->at(0)->beams().at(1)));
#else
    ASSERT_EQ(*(outputTestDataHandle->at(0)->signal_process_vertex()), *(inputTestDataHandle->at(0)->signal_process_vertex()));
    ASSERT_EQ(**(outputTestDataHandle->at(0)->vertices_begin()), **(inputTestDataHandle->at(0)->vertices_begin()));
    ASSERT_EQ(*(outputTestDataHandle->at(0)->beam_particles().first), *(inputTestDataHandle->at(0)->beam_particles().first));
    ASSERT_EQ(*(outputTestDataHandle->at(0)->beam_particles().second), *(inputTestDataHandle->at(0)->beam_particles().second));
#endif
  }

} // <-- namespace SimTesting


int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest( &argc, argv );
  ::testing::AddGlobalTestEnvironment( new SimTesting::GaudiEnvironment );
  return RUN_ALL_TESTS();
}
