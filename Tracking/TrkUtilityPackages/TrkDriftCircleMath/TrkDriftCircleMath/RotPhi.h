/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DCMATH_ROTPHI_H
#define DCMATH_ROTPHI_H

#include <cmath>

#include "TrkDriftCircleMath/LocVec2D.h"

namespace TrkDriftCircleMath {

    class RotPhi {
    public:
        RotPhi(double phi) : m_phi{phi}, m_cosphi{std::cos(phi)}, m_sinphi{std::sin(phi)} {}

        RotPhi(const RotPhi&) = default;
        RotPhi(RotPhi&&) = default;

        RotPhi& operator=(RotPhi&&) = default;
        RotPhi& operator=(const RotPhi&) = default;

        double phi() const { return m_phi; }
        double cosphi() const { return m_cosphi; }
        double sinphi() const { return m_sinphi; }

        double xval(const LocVec2D& lv) const { return cosphi() * lv.x() + sinphi() * lv.y(); }
        double yval(const LocVec2D& lv) const { return -sinphi() * lv.x() + cosphi() * lv.y(); }
        LocVec2D operator*(const LocVec2D& lv) const {return LocVec2D{xval(lv), yval(lv)}; }

        RotPhi inverse() const { return RotPhi(-phi()); }

        void set(double phi) {
            m_phi = phi;
            m_cosphi = std::cos(phi);
            m_sinphi = std::sin(phi);
        }

    private:
        double m_phi{0};
        double m_cosphi{0};
        double m_sinphi{0};
    };

}  // namespace TrkDriftCircleMath
#endif
