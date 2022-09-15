/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "ActsTrkEvent/MultiTrajectory.h"

namespace ActsTrk {
template<>
MultiTrajectory<IsReadWrite>::MultiTrajectory(MultiTrajectory<IsReadWrite>::TrackStateContainerBackendPtr states, 
                                            MultiTrajectory<IsReadWrite>::TrackParametersContainerBackendPtr parameters,
                                            MultiTrajectory<IsReadWrite>::TrackJacobianContainerBackendPtr jacobians, 
                                            MultiTrajectory<IsReadWrite>::TrackMeasurementsContainerBackendPtr measurements )
    : m_trackStates(states),
      m_trackParameters(parameters),
      m_jacobians(jacobians),
      m_measurements(measurements)
{} 

template<>
MultiTrajectory<IsReadOnly>::MultiTrajectory(MultiTrajectory<IsReadWrite>&& rhs) 
    : m_trackStates(rhs.m_trackStates),
      m_trackParameters(rhs.m_trackParameters), 
      m_jacobians(rhs.m_jacobians),
      m_measurements(rhs.m_measurements)
{
    rhs.m_trackStates = nullptr;
    rhs.m_trackParameters = nullptr;
    rhs.m_jacobians = nullptr;
    rhs.m_measurements = nullptr;
} 


/////// MWMW 

    std::size_t
        addTrackState_impl(Acts::TrackStatePropMask mask, std::uint32_t previous) {
        assert(m_trackStates && "Missing Track States backend");
        //auto state = new xAOD::ActsTrackState();
        auto state = new xAOD::TrackState_v1();
        //trackStates().push_back(state);
        // MWMW
        std::cout<<m_trackStates<<std::endl;
        m_trackStates.push_back(state);

        
        //MWMW
        Acts::MultiTrajectoryTraits::IndexType kInvalid = Acts::MultiTrajectoryTraits::kInvalid;
        
        
        
        if (previous >= kInvalid - 1) previous = kInvalid; // fix needed in Acts::MTJ
        trackStates().back()->setPrevious(previous);
        using namespace Acts;

        auto addParam = [this]() -> IndexType {
            trackParameters().push_back(new xAOD::ActsTrackParameters());
            // trackParameters()->back()->resize();
            return trackParameters().size() - 1;
        };

        state->setPredicted(kInvalid);
        if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Predicted)) {
            state->setPredicted(addParam());
        }

        state->setFiltered(kInvalid);
        if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Filtered)) {
            state->setFiltered(addParam());
        }
        state->setSmoothed(kInvalid);
        if (ACTS_CHECK_BIT(mask, TrackStatePropMask::Smoothed)) {
            state->setSmoothed(addParam());
        }
        state->setJacobian(kInvalid);
        state->setUncalibrated(kInvalid);
        state->setCalibrated(kInvalid);
        state->setProjector(kInvalid);

        return m_trackStates->size() - 1;
    }

///////


} // EOF namespace ActsTrk 
