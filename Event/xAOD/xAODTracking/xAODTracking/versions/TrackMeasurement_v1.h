/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
#ifndef XAODTRACKING_VERSIONS_TRACKMEASUREMENT_V1_H
#define XAODTRACKING_VERSIONS_TRACKMEASUREMENT_V1_H
#include <cstdint>
#include "AthLinks/ElementLink.h"
#include "AthContainers/AuxElement.h"
#include "xAODMeasurementBase/UncalibratedMeasurementContainer.h"

#include "EventPrimitives/EventPrimitives.h"


static const SG::AuxElement::Accessor<std::vector<double>> measAcc("meas");
static const SG::AuxElement::Accessor<std::vector<double>> covMatrixAcc("covMatrix");

namespace xAOD {
    /**
     * @brief Track Measurements for Acts MultiTrajectory
     **/

    class TrackMeasurement_v1 : public SG::AuxElement {
    public:
        TrackMeasurement_v1() = default;
        /**
         * access track Measurements vector of const element
         **/
	template<std::size_t measdim = 6>
	  Eigen::Map<const Eigen::Matrix<double, measdim, 1>> measEigen() const {
	  return Eigen::Map<const Eigen::Matrix<double, measdim, 1>>{measAcc(*this).data()}; 
	}
        /**
         * access Measurements of non const element
         **/
	template<std::size_t measdim = 6>
	  Eigen::Map<Eigen::Matrix<double, measdim, 1>> measEigen() {
	  return Eigen::Map<Eigen::Matrix<double, measdim, 1>>{measAcc(*this).data()}; 
	}

        /**
         * access track Measurements as plain vector
         **/
        const std::vector<double>& meas() const;
        /**
         * access set Measurements from plain vector
         **/
        void setMeas( const std::vector<double>& m);


        /**
         * access track covariance matrix (flattened, rows layout) of const element
         **/
	template<std::size_t measdim = 6>
	  Eigen::Map<const Eigen::Matrix<double, measdim, measdim>> covMatrixEigen() const {
	  return Eigen::Map<const Eigen::Matrix<double, measdim, measdim>>{covMatrixAcc(*this).data()};
	}

        /**
         * access track covariance matrix (flattened, rows layout)
         **/
	template<std::size_t measdim = 6>
	  Eigen::Map<Eigen::Matrix<double, measdim, measdim>> covMatrixEigen() {
	  return Eigen::Map<Eigen::Matrix<double, measdim, measdim>> {covMatrixAcc(*this).data()};
	}

        /**
         * access track covariance as plain vector
         **/
        const std::vector<double>& covMatrix() const;
        /**
         * access set covariance from plain vector
         **/
        void setCovMatrix( const std::vector<double>& m);


        /**
         * @brief access the uncalibrated measurement
         * TODO consider bare pointer access
         */
        const ElementLink<UncalibratedMeasurementContainer>& uncalibratedMeasurementLink() const;

        /**
         * @brief return pointer to uncalibrated measurement if the underlying link is valid
         * 
         * @return const UncalibratedMeasurement* or nullptr
         */
        const UncalibratedMeasurement* uncalibratedMeasurement() const;

        /**
         * @brief set uncalibrated measurement
         */
        void setUncalibratedMeasurementLink( const ElementLink<UncalibratedMeasurementContainer>& link);


        // TODO projectors

        /**
         * @brief expands sizes of internal vectors for the data storage
         * ( by default this is 6 for Measurements and 6x6 for Covariance)
         * typically only 5 are used
         */
        void resize(size_t sz = 6);
    };
}
#endif
