/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
#define BOOST_TEST_MODULE MultiTrajectoryBasic_test
#include <boost/test/data/test_case.hpp>
#include <boost/test/included/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
//#include "Acts/EventData/TrackParameters.hpp"

#include "ActsTrkEvent/MultiTrajectory.h"
#include "xAODTracking/TrackStateAuxContainer.h"
#include "xAODTracking/TrackParametersAuxContainer.h"
#include "xAODTracking/TrackJacobianAuxContainer.h"
#include "xAODTracking/TrackMeasurementAuxContainer.h"


namespace {


struct EmptyMTJ { // setup empty MTJ
  EmptyMTJ() {
    trackStateBackend = std::make_unique<xAOD::TrackStateContainer>();
    trackStateBackendAux = std::make_unique<xAOD::TrackStateAuxContainer>();
    trackStateBackend->setStore(trackStateBackendAux.get());

    parametersBackend = std::make_unique<xAOD::TrackParametersContainer>();
    parametersBackendAux = std::make_unique<xAOD::TrackParametersAuxContainer>();
    parametersBackend->setStore(parametersBackendAux.get());

    jacobianBackend = std::make_unique<xAOD::TrackJacobianContainer>();
    jacobianBackendAux = std::make_unique<xAOD::TrackJacobianAuxContainer>();
    jacobianBackend->setStore(jacobianBackendAux.get());

    measurementsBackend = std::make_unique<xAOD::TrackMeasurementContainer>();
    measurementsBackendAux = std::make_unique<xAOD::TrackMeasurementAuxContainer>();
    measurementsBackend->setStore(measurementsBackendAux.get());

    mtj = std::make_unique<ActsTrk::MutableMultiTrajectory>(trackStateBackend.get(), parametersBackend.get(), 
                                                            jacobianBackend.get(), measurementsBackend.get());
    // backends can be shared 
    ro_mtj = std::make_unique<ActsTrk::ConstMultiTrajectory>(trackStateBackend.get(), parametersBackend.get(), 
                                                            jacobianBackend.get(), measurementsBackend.get());

  }
  std::unique_ptr<xAOD::TrackStateContainer> trackStateBackend;
  std::unique_ptr<xAOD::TrackStateAuxContainer> trackStateBackendAux;
  std::unique_ptr<xAOD::TrackParametersContainer> parametersBackend;
  std::unique_ptr<xAOD::TrackParametersAuxContainer> parametersBackendAux;
  std::unique_ptr<xAOD::TrackJacobianContainer> jacobianBackend;
  std::unique_ptr<xAOD::TrackJacobianAuxContainer> jacobianBackendAux;
  std::unique_ptr<xAOD::TrackMeasurementContainer> measurementsBackend;
  std::unique_ptr<xAOD::TrackMeasurementAuxContainer> measurementsBackendAux;


  std::unique_ptr<ActsTrk::MutableMultiTrajectory> mtj;
  std::unique_ptr<ActsTrk::ConstMultiTrajectory> ro_mtj;
};


BOOST_FIXTURE_TEST_CASE(Fill, EmptyMTJ) {
    BOOST_CHECK( mtj->has_backends());
    constexpr auto kMask = Acts::TrackStatePropMask::Predicted;
    auto i0 = mtj->addTrackState(kMask);
    auto i1a = mtj->addTrackState(kMask, i0);
    auto i2a = mtj->addTrackState(kMask, i1a);

    


  std::vector<size_t> act;
  auto collect = [&](auto p) {
    act.push_back(p.index());
    BOOST_CHECK(!p.hasCalibrated());
    BOOST_CHECK(!p.hasFiltered());
    BOOST_CHECK(!p.hasSmoothed());
    BOOST_CHECK(!p.hasJacobian());
    BOOST_CHECK(!p.hasProjector());
  };

  
 
  const std::vector<size_t> exp = { i2a, i1a, i0 };

  
  mtj->visitBackwards(i2a, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());
  // the same test on read only collection TODO, this needs streamlining so we so not repeat identical code
  act.clear();
  ro_mtj->visitBackwards(i2a, collect);
  BOOST_CHECK_EQUAL_COLLECTIONS(act.begin(), act.end(), exp.begin(), exp.end());

}

BOOST_FIXTURE_TEST_CASE(Dynamic_columns, EmptyMTJ) {
    using namespace Acts::HashedStringLiteral;
    BOOST_CHECK( mtj->has_backends());
    mtj->enableDecoration<short>("author");
    constexpr auto kMask = Acts::TrackStatePropMask::Predicted;
    auto i0 = mtj->addTrackState(kMask);
    auto i1 = mtj->addTrackState(kMask, i0);
    auto i2 = mtj->addTrackState(kMask, i1);
    // dynamic column enabled late
    mtj->enableDecoration<float>("mcprob");
    auto ts0 = mtj->getTrackState(i0);
    auto ts1 = mtj->getTrackState(i1);
    auto ts2 = mtj->getTrackState(i2);
    ts0.component<short, "author"_hash>() = 5;
    ts1.component<short, "author"_hash>() = 6;
    ts2.component<short, "author"_hash>() = 4;

    ts0.component<float, "mcprob"_hash>() = 0.5;
    ts1.component<float, "mcprob"_hash>() = 0.9;
    // unset for ts2

    // read them back
    BOOST_CHECK_EQUAL((ts0.component<short,"author"_hash>()), 5);
    BOOST_CHECK_EQUAL((ts1.component<short,"author"_hash>()), 6);
    BOOST_CHECK_EQUAL((ts2.component<short,"author"_hash>()), 4);
    
    BOOST_TEST((ts0.component<float,"mcprob"_hash>()) == 0.5, boost::test_tools::tolerance(0.01));
    BOOST_TEST((ts1.component<float,"mcprob"_hash>()) == 0.9, boost::test_tools::tolerance(0.01));
    BOOST_TEST((ts2.component<float,"mcprob"_hash>()) == 0.0, boost::test_tools::tolerance(0.01));

    BOOST_CHECK_THROW((ts2.component<float,"sth"_hash>()), std::runtime_error);
}
// TODO remaining tests
}