/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef LARHV_EMECHVELECTRODE_H
#define LARHV_EMECHVELECTRODE_H

class EMECHVModule;

class EMECHVElectrode
{
 public:
  // The ownership of the EMBHVModule pointer NOT transferred to EMECHVElectrode
  EMECHVElectrode(const EMECHVModule* module, unsigned int iElectrode);
  ~EMECHVElectrode();

  unsigned int getElectrodeIndex() const;

  // returns a pointer to the module that owns this electrode.
  const EMECHVModule& getModule() const;

  double getPhi() const;

  // HV Status
  bool hvOn(int iGap) const;
    
  // Voltage
  double voltage(int iGap) const;
  
  // Current
  double current(int iGap) const;
  
  // Voltage and current at the same time
  void voltage_current(int iGap, double& v, double& i) const;

  // HVLine no
  int hvLineNo(int iGap) const;

 private:
  EMECHVElectrode(const EMECHVElectrode& right);
  EMECHVElectrode& operator=(const EMECHVElectrode& right);

  class Clockwork;
  Clockwork *m_c;
};

#endif
