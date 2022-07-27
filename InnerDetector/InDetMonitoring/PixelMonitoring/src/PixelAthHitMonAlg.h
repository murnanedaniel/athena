/*
   Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
 */

#ifndef PIXELATHHITMONTOOL_H
#define PIXELATHHITMONTOOL_H

#include "PixelAthMonitoringBase.h"
#include "InDetConditionsSummaryService/IInDetConditionsTool.h"
#include "InDetRawData/PixelRDO_Container.h"
#include "PixelReadoutGeometry/IPixelReadoutManager.h"

class PixelID;

class PixelAthHitMonAlg: public PixelAthMonitoringBase {
public:
  PixelAthHitMonAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~PixelAthHitMonAlg() = default;
  virtual StatusCode initialize() override;
  virtual StatusCode fillHistograms(const EventContext& ctx) const override;
  std::string findComponentString(int bec, int ld) const;
private:

  SG::ReadHandleKey<PixelRDO_Container> m_pixelRDOName {
    this, "RDOName", "PixelRDOs", "rdo data key"
  };

  bool m_doOnline {};
  bool m_doLumiBlock {};
  bool m_doLowOccupancy {};
  bool m_doHighOccupancy {};
  bool m_doHeavyIonMon {};
  bool m_doFEPlots {};
};
#endif
