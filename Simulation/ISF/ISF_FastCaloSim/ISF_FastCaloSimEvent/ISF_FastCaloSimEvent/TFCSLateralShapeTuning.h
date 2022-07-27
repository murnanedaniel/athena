/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSLateralShapeTuning_h
#define TFCSLateralShapeTuning_h

/*Parent class*/
#include "ISF_FastCaloSimEvent/TFCSLateralShapeParametrizationHitBase.h"
/*Truth state include*/
#include "ISF_FastCaloSimEvent/TFCSTruthState.h"
/*Particle data includes*/
#include "HepPDT/ParticleData.hh"
#include "HepPDT/ParticleDataTable.hh"
/*TFile include*/
#include <TFile.h>
/*Pieciewise linear interpolation class include*/
#include "ISF_FastCaloSimEvent/TFCSEnergyInterpolationPiecewiseLinear.h"


class TFCSLateralShapeTuning: public TFCSLateralShapeParametrizationHitBase {

  typedef std::map<std::string, TFCSEnergyInterpolationPiecewiseLinear*> interpolationMap;

  public:

    /// Constructor
    TFCSLateralShapeTuning(const char* name=nullptr, const char* title=nullptr);
    // Destructor
    ~TFCSLateralShapeTuning();
    // Simulation call
    virtual FCSReturnCode simulate_hit(Hit& hit, TFCSSimulationState& simulstate, const TFCSTruthState* truth, const TFCSExtrapolationState* extrapol) override;
    // Initialization from model parameter file
    FCSReturnCode initFromModelFile(const std::string& pathToModelParameters, int intMinEta, int intMaxEta);
    // Initialization from interpolation map
    FCSReturnCode initFromMap(const interpolationMap&);
    // Helper function for strip layer model
    double getSeriesScalingFactor(double a0, double a1, double a2, double a3, double distToShowerCenter);

  private:

    // mapping between the parameter names of the model and its pieciewise linear interpolation obect
    interpolationMap m_parameterInterpol;

    ClassDefOverride(TFCSLateralShapeTuning, 1)  //TFCSLateralShapeTuning
};

#endif
