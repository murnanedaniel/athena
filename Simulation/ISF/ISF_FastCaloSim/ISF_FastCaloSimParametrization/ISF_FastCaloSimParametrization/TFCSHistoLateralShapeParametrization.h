/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHistoLateralShapeParametrization_h
#define TFCSHistoLateralShapeParametrization_h

#include "ISF_FastCaloSimParametrization/TFCSLateralShapeParametrization.h"

#include "TH2.h"
#include "TRandom3.h"


class TFCSHistoLateralShapeParametrization:public TFCSLateralShapeParametrization {
public:
  TFCSHistoLateralShapeParametrization(const char* name=0, const char* title=0);

  // simulated one hit position with weight that should be put into simulstate
  // sometime later all hit weights should be resacled such that their final sum is simulstate->E(sample)
  // someone also needs to map all hits into cells
  virtual void simulate(TFCSSimulationState& simulstate,const TFCSTruthState* truth, const TFCSExtrapolationState* extrapol);

  // Init and fill sigma
  bool Initialize(TH2* hist);
  bool Initialize(const char* filepath, const char* histname);

private:
  // simple shape information should be stored as private member variables here

  TH2* m_hist;

  TRandom3 *m_rnd;

  ClassDef(TFCSHistoLateralShapeParametrization,1)  //TFCSHistoLateralShapeParametrization
};

#if defined(__MAKECINT__)
#pragma link C++ class TFCSHistoLateralShapeParametrization+;
#endif

#endif
