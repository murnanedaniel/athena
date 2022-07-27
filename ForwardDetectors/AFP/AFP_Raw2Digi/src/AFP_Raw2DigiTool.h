/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#ifndef AFP_RAW2DIGITOOL_H
#define AFP_RAW2DIGITOOL_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "AFP_Raw2Digi/IAFP_Raw2DigiTool.h"
#include "AFP_RawEv/AFP_RawContainer.h"
#include "xAODForward/AFPToFHit.h"
#include "xAODForward/AFPToFHitContainer.h"
#include "xAODForward/AFPSiHitContainer.h"
#include <string>
#include "TF1.h"

class AFP_Raw2DigiTool : public extends<AthAlgTool, IAFP_Raw2DigiTool>
{
public:
  AFP_Raw2DigiTool(const std::string &type, const std::string &name, const IInterface *parent);

  /// Does nothing 
  virtual ~AFP_Raw2DigiTool() override;

  /// Does nothing
  virtual StatusCode initialize() override;

  /// Creates xAOD for silicon detector
  StatusCode recoSiHits(const EventContext &ctx) const override;

  /// Creates xAOD for time-of-flight detector
  StatusCode recoToFHits(const EventContext &ctx) const override;

  StatusCode recoAll(const EventContext &ctx) const override;

  
  /// Does nothing
  virtual StatusCode finalize() override;

  /// Method that decodes raw information about time-over-threshold to
  /// number of clock ticks
  unsigned int decodeTimeOverThresholdSi (const unsigned int input, const unsigned int discConfig) const;
  
protected:
  SG::ReadHandleKey<AFP_RawContainer> m_rawDataContainerName{this, "rawDataContainerName", "AFP_RawData"};
  SG::WriteHandleKey<xAOD::AFPSiHitContainer> m_AFPSiHitsContainerName{this, "AFPSiHitsContainerName", "AFPSiHitContainer"};
  SG::WriteHandleKey<xAOD::AFPToFHitContainer> m_AFPHitsContainerNameToF{this, "AFPHitsContainerNameToF", "AFPToFHitContainer"};

  /// @brief Factor converting signal to time
  ///
  /// The value of the factor is 25/1024 nanoseconds
  static constexpr double s_timeConversionFactor = 25./1024.;

  /// @brief Factor converting pulse length to time
  ///
  /// The value of the factor is 0.521 nanoseconds
  static constexpr double s_pulseLengthFactor = 0.521;

  /// @brief Function that transforms time-over-threshold to charge
  ///
  /// Transformation function can be set in steering cards
  Gaudi::Property<std::string> m_totToChargeTransfExpr{this, "TotToChargeTransfExpr", "1909 + x*363 + x*x*141", "Function that transforms time-over-threshold to charge"};
  Gaudi::Property<std::string> m_totToChargeTransfName{this, "TotToChargeTransfName", "TotToChargeTransfFunction", "Name of the function that transforms time-over-threshold to charge"};
  TF1 m_totToChargeTransformation;	

  /// Method that creates a new AFPToFHit and sets it valus according to #data
  void newXAODHitToF (xAOD::AFPToFHitContainer* tofHitContainer, const AFP_ToFRawCollection& collection, const AFP_ToFRawData& data, const EventContext& ctx) const;
  
  /// Method that creates a new AFPSiHit and sets it valus according to #data
  void newXAODHitSi (xAOD::AFPSiHitContainer* xAODSiHit, const AFP_SiRawCollection& collection, const AFP_SiRawData& data) const;

  /// @brief Method mapping hptdcID and hptdcChannel to train ID and bar in train ID
  ///
  /// The method requires that hptdcID and hptdcChannel are set in the
  /// tofHit passed as argument.  Mapping is implemented according to
  /// https://twiki.cern.ch/twiki/bin/view/Atlas/AFPHPTDC#Channel_Mapping
  void setBarAndTrainID(xAOD::AFPToFHit* tofHit, const EventContext& ctx) const;
};
#endif 
