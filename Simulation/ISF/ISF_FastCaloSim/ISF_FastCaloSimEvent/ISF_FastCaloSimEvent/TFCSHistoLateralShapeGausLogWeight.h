/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSHistoLateralShapeGausLogWeight_h
#define TFCSHistoLateralShapeGausLogWeight_h

#include "ISF_FastCaloSimEvent/TFCSHistoLateralShapeWeight.h"
#ifdef USE_GPU
#include "ISF_FastCaloGpu/LoadGpuHist.h"
#endif

class TH1;

class TFCSHistoLateralShapeGausLogWeight:public TFCSHistoLateralShapeWeight {
public:
  TFCSHistoLateralShapeGausLogWeight(const char* name=nullptr, const char* title=nullptr);
  virtual ~TFCSHistoLateralShapeGausLogWeight();

  /// weight the energy of one hit in order to generate fluctuations
  virtual FCSReturnCode simulate_hit(Hit& hit,TFCSSimulationState& simulstate,const TFCSTruthState* truth, const TFCSExtrapolationState* extrapol) override;

#ifdef USE_GPU
  void    set_d_Hist( FH1D* hf_ptr ) { m_d_Hist = hf_ptr; };
  const FH1D*       d_Hist() { return m_d_Hist; };
  void             LoadHist();
  LoadGpuHist* LdFH() { return m_LdFH; };
#endif

protected:

  ClassDefOverride(TFCSHistoLateralShapeGausLogWeight,1)  //TFCSHistoLateralShapeGausLogWeight

private:

#ifdef USE_GPU
  FH1D*             m_d_Hist = nullptr;
  LoadGpuHist* m_LdFH        = nullptr;
#endif
};

#endif
