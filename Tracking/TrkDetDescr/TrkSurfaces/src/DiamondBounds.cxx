/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// DiamondBounds.cxx, (c) ATLAS Detector Software
///////////////////////////////////////////////////////////////////

// Trk
#include "TrkSurfaces/DiamondBounds.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
// STD
#include <cmath>
#include <iomanip>
#include <iostream>

// default constructor
Trk::DiamondBounds::DiamondBounds()
  : m_boundValues(DiamondBounds::bv_length, 0.)
  , m_alpha1(0.)
  , m_alpha2(0.)
{}

// constructor from arguments I
Trk::DiamondBounds::DiamondBounds(double minhalex, double medhalex, double maxhalex, double haley1, double haley2)
  : m_boundValues(DiamondBounds::bv_length, 0.)
  , m_alpha1(0.)
  , m_alpha2(0.)
{
  m_boundValues[DiamondBounds::bv_minHalfX] = minhalex;
  m_boundValues[DiamondBounds::bv_medHalfX] = medhalex;
  m_boundValues[DiamondBounds::bv_maxHalfX] = maxhalex;
  m_boundValues[DiamondBounds::bv_halfY1] = haley1;
  m_boundValues[DiamondBounds::bv_halfY2] = haley2;
  if (minhalex > maxhalex)
    swap(m_boundValues[DiamondBounds::bv_minHalfX], m_boundValues[DiamondBounds::bv_maxHalfX]);
  DiamondBounds::initCache();
}

bool
Trk::DiamondBounds::operator==(const Trk::SurfaceBounds& sbo) const
{
  // check the type first not to compare apples with oranges
  const Trk::DiamondBounds* diabo = dynamic_cast<const Trk::DiamondBounds*>(&sbo);
  if (!diabo)
    return false;
  return (m_boundValues == diabo->m_boundValues);
}

// checking if inside bounds (Full symmetrical Diamond)
bool
Trk::DiamondBounds::inside(const Amg::Vector2D& locpo, double tol1, double tol2) const
{
  return this->insideFull(locpo, tol1, tol2);
}

// checking if inside bounds (Full symmetrical Diamond)
bool
Trk::DiamondBounds::insideFull(const Amg::Vector2D& locpo, double tol1, double tol2) const
{
  // validity check
  if (!m_boundValues[DiamondBounds::bv_halfY1] && !m_boundValues[DiamondBounds::bv_minHalfX]) return false;

  // quick False (radial direction)
  if (locpo[Trk::locY] < -2. * m_boundValues[DiamondBounds::bv_halfY1] - tol2) return false;
  if (locpo[Trk::locY] >  2. * m_boundValues[DiamondBounds::bv_halfY2] + tol2) return false;

  double absX = std::abs(locpo[Trk::locX]);

  // quick False (transverse directon)
  if (absX > m_boundValues[DiamondBounds::bv_medHalfX] + tol1) return false;

  // quick True
  if (absX < std::min(m_boundValues[DiamondBounds::bv_minHalfX], m_boundValues[DiamondBounds::bv_maxHalfX]) + tol1) return true;

  /** use basic calculation of a straight line */
  const double& halfBaseUp = locpo[Trk::locY] < 0 ? m_boundValues[DiamondBounds::bv_medHalfX] : m_boundValues[DiamondBounds::bv_maxHalfX];
  const double& halfBaseLo = locpo[Trk::locY] < 0 ? m_boundValues[DiamondBounds::bv_minHalfX] : m_boundValues[DiamondBounds::bv_medHalfX];
  const double& halfH      = locpo[Trk::locY] < 0 ? m_boundValues[DiamondBounds::bv_halfY1]   : m_boundValues[DiamondBounds::bv_halfY2];
  double k = halfH ? 0.5*(halfBaseUp - halfBaseLo)/halfH : 0.;
  double sign = (k < 0) ? -1. : 1.;   
  return (absX - tol1 <= m_boundValues[DiamondBounds::bv_medHalfX] + k * (locpo[Trk::locY] + sign*tol2));

}

// opening angle in point A
double
Trk::DiamondBounds::alpha1() const
{
  return m_alpha1;
}

// opening angle in point A'
double
Trk::DiamondBounds::alpha2() const
{
  return m_alpha2;
}

double
Trk::DiamondBounds::minDistance(const Amg::Vector2D& pos) const
{
  const int Np = 6;

  double y1 = 2. * m_boundValues[DiamondBounds::bv_halfY1];
  double y2 = 2. * m_boundValues[DiamondBounds::bv_halfY2];

  double X[6] = { -m_boundValues[DiamondBounds::bv_minHalfX], -m_boundValues[DiamondBounds::bv_medHalfX],
                  -m_boundValues[DiamondBounds::bv_maxHalfX], m_boundValues[DiamondBounds::bv_maxHalfX],
                  m_boundValues[DiamondBounds::bv_medHalfX],  m_boundValues[DiamondBounds::bv_minHalfX] };
  double Y[6] = { -y1, 0., y2, y2, 0., -y1 };

  double dm = 1.e+20;
  double Ao = 0.;
  bool in = true;

  for (int i = 0; i != Np; ++i) {

    int j = (i == Np-1 ? 0 : i+1);

    double x = X[i] - pos[0];
    double y = Y[i] - pos[1];
    double dx = X[j] - X[i];
    double dy = Y[j] - Y[i];
    double A = x * dy - y * dx;
    double S = -(x * dx + y * dy);

    if (S <= 0.) {
      double d = x * x + y * y;
      if (d < dm)
        dm = d;
    } else {
      double a = dx * dx + dy * dy;
      if (S <= a) {
        double d = (A * A) / a;
        if (d < dm)
          dm = d;
      }
    }
    if (i && in && Ao * A < 0.)
      in = false;
    Ao = A;
  }
  if (in)
    return -sqrt(dm);
  return sqrt(dm);
}

// ostream operator overload

MsgStream&
Trk::DiamondBounds::dump(MsgStream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Trk::DiamondBounds:  (minHlenghtX, medHlengthX, maxHlengthX, hlengthY1, hlengthY2 ) = ";
  sl << "(" << m_boundValues[DiamondBounds::bv_minHalfX] << ", " << m_boundValues[DiamondBounds::bv_medHalfX] << ", "
     << m_boundValues[DiamondBounds::bv_maxHalfX] << ", " << m_boundValues[DiamondBounds::bv_halfY1] << ", "
     << m_boundValues[DiamondBounds::bv_halfY2] << ")";
  sl << std::setprecision(-1);
  return sl;
}

std::ostream&
Trk::DiamondBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Trk::DiamondBounds:  (minHlenghtX, medHlengthX, maxHlengthX, hlengthY1, hlengthY2 ) = ";
  sl << "(" << m_boundValues[DiamondBounds::bv_minHalfX] << ", " << m_boundValues[DiamondBounds::bv_medHalfX] << ", "
     << m_boundValues[DiamondBounds::bv_maxHalfX] << ", " << m_boundValues[DiamondBounds::bv_halfY1] << ", "
     << m_boundValues[DiamondBounds::bv_halfY2] << ")";
  sl << std::setprecision(-1);
  return sl;
}
