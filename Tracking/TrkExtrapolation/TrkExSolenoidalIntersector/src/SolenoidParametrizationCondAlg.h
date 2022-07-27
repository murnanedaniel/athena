/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/ 

#ifndef _SolenoidParametrizationCondAlg_H_
#define _SolenoidParametrizationCondAlg_H_

#include "AthenaBaseComps/AthReentrantAlgorithm.h"

#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"

#include "TrkExSolenoidalIntersector/SolenoidParametrization.h"
#include "MagFieldConditions/AtlasFieldCacheCondObj.h"

#include "GaudiKernel/ICondSvc.h"

namespace Trk
{

class SolenoidParametrizationCondAlg : public AthReentrantAlgorithm
{
 public:
  SolenoidParametrizationCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~SolenoidParametrizationCondAlg() = default;
  virtual StatusCode initialize() override final;
  virtual StatusCode execute(const EventContext &ctx) const override final;
  virtual StatusCode finalize() override final;
  virtual bool isReEntrant() const override final { return false; }

 private:
   SG::ReadCondHandleKey<AtlasFieldCacheCondObj> m_fieldCacheCondObjInputKey
      {this, "AtlasFieldCacheCondObj", "fieldCondObj", "Name of the Magnetic Field conditions object key"};

   SG::WriteCondHandleKey<SolenoidParametrization> m_writeKey
      { this, "WriteKey", "SolenoidParametrization"};

   ServiceHandle<ICondSvc> m_condSvc
      { this, "CondSvc", "CondSvc"};
};
}
#endif // SCT_CONDITIONSPARAMETERCONDALG
