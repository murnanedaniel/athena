/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef GENERATORFILTERSVBFMJJINTERVALFILTER_H
#define GENERATORFILTERSVBFMJJINTERVALFILTER_H

#include "GeneratorModules/GenFilter.h"
#include "JetEvent/JetCollection.h"
#include "JetEvent/Jet.h"

class IAtRndmGenSvc;


class VBFMjjIntervalFilter : public GenFilter {
public:

  VBFMjjIntervalFilter(const std::string& name, ISvcLocator* pSvcLocator);
  virtual StatusCode filterInitialize();
  virtual StatusCode filterEvent();

private:

  double m_olapPt;
  double m_yMax;                           // Rapidity acceptance
  double m_pTavgMin;                       // Required average dijet pT
  std::string m_TruthJetContainerName;     // Name of the truth jet container

  ServiceHandle<IAtRndmGenSvc> m_rand;     // Random number generator

  long m_total;                            // Total number of events tested
  long m_passed;                           // Number of events passing all cuts
  //long m_outsideAcceptance;                // Number of events failing rapidity acceptance cuts

  double m_norm;                           // Normalization for weights
  //double m_high;                           // High-side function level
  //bool m_doShape;                          // Attempt to flatten the dY distribution
  double m_prob0;
  double m_prob1;
  double m_prob2low;
  double m_prob2high;
  double m_mjjlow;
  double m_mjjhigh;
  double m_alpha;

  bool checkOverlap(double, double, std::vector<HepMC::GenParticle*>);


public:

  double getEventWeight(JetCollection *jets);
};

#endif
