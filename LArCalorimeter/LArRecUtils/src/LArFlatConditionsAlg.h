//Dear emacs, this is -*- C++ -*- 

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef LARFLATCONDITIONSALG_H
#define LARFLATCONDITIONSALG_H

#include "AthenaBaseComps/AthAlgorithm.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"

#include "AthenaPoolUtilities/CondAttrListCollection.h"

template<class T>
class LArFlatConditionsAlg: public AthAlgorithm {
 public:
  using AthAlgorithm::AthAlgorithm;
  ~LArFlatConditionsAlg() = default;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

 private:
  SG::ReadCondHandleKey<CondAttrListCollection> m_readKey{this,"ReadKey","","Key of the input CDO (AttrListCollection)"};
  SG::WriteCondHandleKey<T>  m_writeKey{this,"WriteKey","","Key of the output LArXYZFlat CDO"};
};

#include "LArFlatConditionsAlg.icc"




#include "LArCOOLConditions/LArPedestalFlat.h"
typedef LArFlatConditionsAlg<LArPedestalFlat> LArCondAlgPedestalFlat;

#include "LArCOOLConditions/LArAutoCorrSC.h"
typedef LArFlatConditionsAlg<LArAutoCorrSC> LArCondAlgAutoCorrSC;

#include "LArCOOLConditions/LArDAC2uAFlat.h"
typedef LArFlatConditionsAlg<LArDAC2uAFlat> LArCondAlgDAC2uAFlat;

#include "LArCOOLConditions/LArDAC2uASC.h"
typedef LArFlatConditionsAlg<LArDAC2uASC> LArCondAlgDAC2uASC;

#include "LArCOOLConditions/LArHVScaleCorrFlat.h"
typedef LArFlatConditionsAlg<LArHVScaleCorrFlat> LArCondAlgHVScaleCorrFlat;

#include "LArCOOLConditions/LArHVScaleCorrSC.h"
typedef LArFlatConditionsAlg<LArHVScaleCorrSC> LArCondAlgHVScaleCorrSC;

#include "LArCOOLConditions/LArMinBiasSC.h"
typedef LArFlatConditionsAlg<LArMinBiasSC> LArCondAlgMinBiasSC;

#include "LArCOOLConditions/LArMinBiasAverageSC.h"
typedef LArFlatConditionsAlg<LArMinBiasAverageSC> LArCondAlgMinBiasAverageSC;

#include "LArCOOLConditions/LArMphysOverMcalFlat.h"
typedef LArFlatConditionsAlg<LArMphysOverMcalFlat> LArCondAlgMphysOverMcalFlat;

#include "LArCOOLConditions/LArMphysOverMcalSC.h"
typedef LArFlatConditionsAlg<LArMphysOverMcalSC> LArCondAlgMphysOverMcalSC;

#include "LArCOOLConditions/LArNoiseSC.h"
typedef LArFlatConditionsAlg<LArNoiseSC> LArCondAlgNoiseSC;

#include "LArCOOLConditions/LArOFCFlat.h"
typedef LArFlatConditionsAlg<LArOFCFlat> LArCondAlgOFCFlat;

#include "LArCOOLConditions/LArPedestalSC.h"
typedef LArFlatConditionsAlg<LArPedestalSC> LArCondAlgPedestalSC;

#include "LArCOOLConditions/LArRampFlat.h"
typedef LArFlatConditionsAlg<LArRampFlat> LArCondAlgRampFlat;

#include "LArCOOLConditions/LArRampSC.h"
typedef LArFlatConditionsAlg<LArRampSC> LArCondAlgRampSC;

#include "LArCOOLConditions/LArShapeFlat.h"
typedef LArFlatConditionsAlg<LArShapeFlat> LArCondAlgShapeFlat;

#include "LArCOOLConditions/LArShapeSC.h"
typedef LArFlatConditionsAlg<LArShapeSC> LArCondAlgShapeSC;

#include "LArCOOLConditions/LArfSamplSC.h"
typedef LArFlatConditionsAlg<LArfSamplSC> LArCondAlgfSamplSC;

#include "LArCOOLConditions/LAruA2MeVFlat.h"
typedef LArFlatConditionsAlg<LAruA2MeVFlat> LArCondAlguA2MeVFlat;

#include "LArCOOLConditions/LAruA2MeVSC.h"
typedef LArFlatConditionsAlg<LAruA2MeVSC> LArCondAlguA2MeVSC;


#endif
