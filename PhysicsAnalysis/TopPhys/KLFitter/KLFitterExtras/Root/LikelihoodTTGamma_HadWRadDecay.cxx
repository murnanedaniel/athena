/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "KLFitterExtras/LikelihoodTTGamma_HadWRadDecay.h" 
#include <vector>

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma_HadWRadDecay::LikelihoodTTGamma_HadWRadDecay()
{
  // calls base class constructor
}

// --------------------------------------------------------- 
KLFitter::LikelihoodTTGamma_HadWRadDecay::~LikelihoodTTGamma_HadWRadDecay()
{
  // calls base class destructor
}

// --------------------------------------------------------- 
int KLFitter::LikelihoodTTGamma_HadWRadDecay::CalculateLorentzVectors(std::vector <double> const& parameters)
{
  KLFitter::LikelihoodTTGamma::CalculateLorentzVectors(parameters);

  // the hadronic W
  *(fParticlesModel->Boson(0)) += *(fParticlesModel->Photon(0));

  // the hadronic top
  *(fParticlesModel->Parton(4)) += *(fParticlesModel->Photon(0));

  // no error 
  return 1; 
}
