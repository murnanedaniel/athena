/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef XAODTRACKING_VERSIONS_TRACKPARAMETERS_V1_H
#define XAODTRACKING_VERSIONS_TRACKPARAMETERS_V1_H
#include <cstdint>
#include "AthContainers/AuxElement.h"
#include "EventPrimitives/EventPrimitives.h"

namespace xAOD {
    /**
     * @brief Track Parameters for Acts MultiTrajectory
     **/

    class TrackParameters_v1 : public SG::AuxElement {
    public:
        // plain std::vectors are mapped to Eigen vectors or matrices
        // TODO a common location for this types
        using ConstVectorMap = Eigen::Map<const Eigen::Matrix<double, 6, 1>>;
        using VectorMap = Eigen::Map<Eigen::Matrix<double, 6, 1>>;
        using ConstMatrixMap = Eigen::Map<const Eigen::Matrix<double, 6, 6>>;
        using MatrixMap = Eigen::Map<Eigen::Matrix<double, 6, 6>>;

        TrackParameters_v1() = default;

        /**
         * access track parameters vector of const element
         **/
        ConstVectorMap parameters() const;
        /**
         * access parameters of non const element
         **/
        VectorMap parameters();

        /**
         * access track covariance matrix (flattened, rows layout) of const element
         **/
        ConstMatrixMap covariance() const;

        /**
         * access track covariance matrix (flattened, rows layout)
         **/
        MatrixMap covariance();

        /**
         * @brief expands sizes of internal vectors for the data storage
         * ( by default this is 6 for parameters and 6x6 for covariance)
         */
        void resize(size_t sz = 6);
    };
}
#endif