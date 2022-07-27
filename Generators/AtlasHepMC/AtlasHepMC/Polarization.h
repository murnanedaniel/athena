/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
/* Author: Andrii Verbytskyi andrii.verbytskyi@mpp.mpg.de */

#ifndef ATLASHEPMC_POLARIZATION_H
#define ATLASHEPMC_POLARIZATION_H
#ifdef HEPMC3
#include "HepMC3/GenParticle.h"
#include "HepMC3/Attribute.h"
namespace HepMC {
class Polarization {
public:
    Polarization( const double a=0.0, const double b=0.0,
                  bool is_defined = false)
      : m_theta(a), m_phi(b), m_is_defined(is_defined) {};
    ~Polarization() {};
    double theta() const  { return m_theta;}
    double phi() const { return m_phi;}
    bool is_defined() const { return m_is_defined; }
    inline bool operator ==(const Polarization &a) { return (std::abs(a.phi() - phi()) <= 0.00001) && (std::abs(a.theta() - theta()) <= 0.00001);}
private:
    double m_theta;
    double m_phi;
    bool m_is_defined;
};
inline Polarization polarization(HepMC3::GenParticlePtr a) {
    std::shared_ptr<HepMC3::DoubleAttribute> phi_A =a->attribute<HepMC3::DoubleAttribute>("phi");
    std::shared_ptr<HepMC3::DoubleAttribute> theta_A=a->attribute<HepMC3::DoubleAttribute>("theta");
    double phi=(phi_A?phi_A->value():0.0);
    double theta=(theta_A?theta_A->value():0.0);
    bool is_defined = phi_A && theta_A;
    return Polarization(theta,phi,is_defined);
}
inline Polarization polarization(HepMC3::ConstGenParticlePtr a) {
    std::shared_ptr<HepMC3::DoubleAttribute> phi_A =a->attribute<HepMC3::DoubleAttribute>("phi");
    std::shared_ptr<HepMC3::DoubleAttribute> theta_A=a->attribute<HepMC3::DoubleAttribute>("theta");
    double phi=(phi_A?phi_A->value():0.0);
    double theta=(theta_A?theta_A->value():0.0);
    bool is_defined = phi_A && theta_A;
    return Polarization(theta,phi,is_defined);
}
template<class T> void  set_polarization( T a,  Polarization b) {
    a->add_attribute("phi",std::make_shared<HepMC3::DoubleAttribute>(b.phi()));
    a->add_attribute("theta",std::make_shared<HepMC3::DoubleAttribute>(b.theta()));
}
}
#else
#include "HepMC/Polarization.h"
namespace HepMC {
template<class T> void set_polarization( T a,  Polarization b) {
    a->set_polarization(b);
}
template<class T>  Polarization polarization(T a) {
    return a->polarization();
}
}
#endif
#endif
