/*
   Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
 */

#ifndef TRIGAFPHYPO_TRIGAFPDIJETCOMBOHYPOTOOL_H
#define TRIGAFPHYPO_TRIGAFPDIJETCOMBOHYPOTOOL_H

#include <string>

#include "IAFPProtonTransportTool.h"
#include "AFPProtonTransportTool.h"

#include "DecisionHandling/ComboHypoToolBase.h"
#include "Gaudi/Property.h"

#include "xAODForward/AFPTrackContainer.h"

#include "AthenaMonitoringKernel/Monitored.h"
#include "AthenaMonitoringKernel/GenericMonitoringTool.h"

class TrigAFPDijetComboHypoTool : public ComboHypoToolBase {
 public:
  
  TrigAFPDijetComboHypoTool(const std::string& type, const std::string& name, const IInterface* parent);
  
  virtual ~TrigAFPDijetComboHypoTool();

  virtual StatusCode initialize() override;

private:

  virtual bool executeAlg(const std::vector<Combo::LegDecision>& combination) const override;

  ///Boolean corresponding to the decision to record the event based on the selection being met.
  //HLT::Identifier m_decisionId;
 
  //Proton transport parameterizaton file for A side
  Gaudi::Property<std::string> m_protonTransportParamFileNameA{this,"protonTransportParamFileName1","final_parameterization_b1.txt","Name of the file with proton transport parameterization for A side"};
  //Proton transport parameterizaton file for C side 
  Gaudi::Property<std::string> m_protonTransportParamFileNameC{this,"protonTransportParamFileName2","final_parameterization_b2.txt","Name of the file with proton transport parameterization C side"};

  //Proton transport parameterizaton object used in the proton position prediction for ATLAS A side
  ToolHandle<IAFPProtonTransportTool> m_transportBeamA;
  ///Proton transport parameterizaton object used in the proton position prediction for ATLAS C side
  ToolHandle<IAFPProtonTransportTool> m_transportBeamC;

  ///Threshold for the radius distance between parameterization and measurements 
  Gaudi::Property<float> m_maxProtonDist{this,"maxProtonDist",2.0,"Cut on distance between measured and extrapolated proton"};
  ///Threshold for the x distance between parameterization and measurements 
  Gaudi::Property<float> m_maxProtonDiff_x{this,"maxProtonDiff_x",2.5,"Cut on x difference between measured and extrapolated proton"};
  ///Threshold for the y distance between parameterization and measurements 
  Gaudi::Property<float> m_maxProtonDiff_y{this,"maxProtonDiff_y",2.5,"Cut on y difference between measured and extrapolated proton"};

  Gaudi::Property<float> m_beamEnergy{this,"beamEnergy",6500.0,"Energy of one beam i.e. half of centre-of-mass energy"};       ///< energy of one beam i.e. half of centre-of-mass energy
  float m_totalEnergy = 2*m_beamEnergy;      ///< beams centre-of-mass energy 2*#m_beamEnergy

  ///Shift in x position between parameterization and measurements
  Gaudi::Property<float> m_protonPosShift_x{this,"protonPosShift_x",0.0,"Expected difference in x position between measured and predicted proton"};
  ///Shif in y position between parameterization and measurements
  Gaudi::Property<float> m_protonPosShift_y{this,"protonPosShift_y",0.0,"Expected difference in y position between measured and predicted proton"};

  ///Beam alignment corrections in the x position of ATLAS A side
  Gaudi::Property<float> m_alignmentCorrection_nearA{this,"alignmentCorrection_nearA",0.0,"Shift due to alignment in x for protons in near station A"};
  ///Beam alignment corrections in the x position of ATLAS C side
  Gaudi::Property<float> m_alignmentCorrection_nearC{this,"alignmentCorrection_nearC",0.0,"Shift due to alignment in x for protons in near station C"};

  ///Variable to convert from MeV to GeV
  const float m_GeV = 0.001;

  SG::ReadHandleKey< xAOD::AFPTrackContainer > m_AFPtrackCollectionReadKey {this, "AFPTrackContainer", "HLT_AFPTrackContainer", "xAOD AFP track collection"};
  
  // monitoring
  ToolHandle<GenericMonitoringTool> m_monTool { this, "MonTool", "", "Monitoring tool" };

};

#endif
