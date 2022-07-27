// -*- C++ -*-
/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef HGTD_CONDITIONSALGORITHMS_HGTD_DETECTORELEMENTCONDALG_H
#define HGTD_CONDITIONSALGORITHMS_HGTD_DETECTORELEMENTCONDALG_H

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "GeoPrimitives/GeoPrimitives.h"
#include "HGTD_ReadoutGeometry/HGTD_DetectorElementCollection.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "StoreGate/CondHandleKeyArray.h"


class HGTD_DetectorManager;

class HGTD_DetectorElementCondAlg : public AthReentrantAlgorithm
{
 public:
  HGTD_DetectorElementCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~HGTD_DetectorElementCondAlg() override = default;

  virtual StatusCode initialize() override final;
  virtual StatusCode execute(const EventContext& ctx) const override final;
  virtual bool isReEntrant() const override final { return false; }

 private:
  SG::WriteCondHandleKey<InDetDD::HGTD_DetectorElementCollection> m_writeKey
  {this, "WriteKey", "HGTD_DetectorElementCollection", "Key of output HGTD_DetectorElementCollection for HGTD"};

  StringProperty m_detManagerName{this, "DetManagerName", "HGTD", "Name of the DeterctorManager to retrieve"};
  const HGTD_DetectorManager* m_detManager{nullptr};
};

#endif // HGTD_CONDITIONSALGORITHMS_HGTD_DETECTORELEMENTCONDALG_H
