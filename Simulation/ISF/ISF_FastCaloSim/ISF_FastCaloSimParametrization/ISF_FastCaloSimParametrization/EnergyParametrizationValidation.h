/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TFCSSimulationState_h
#define EnergyParametrizationValidation_h

class EnergyParametrizationValidation
{
  public:
    //EnergyParametrizationValidation();
    //virtual ~EnergyParametrizationValidation(){};
    static void autozoom(TH1D* h1, TH1D* h2, double &min, double &max, double &rmin, double &rmax);
    static TH1D* refill(TH1D* h_in,double min, double max, double rmin, double rmax);

  private:
  ClassDef(EnergyParametrizationValidation,1)
};

#endif
