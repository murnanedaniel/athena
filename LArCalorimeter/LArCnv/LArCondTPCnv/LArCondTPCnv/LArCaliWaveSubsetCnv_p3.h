/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef LARCONDTPCNV_LARcaliWAVEsubsetCNV_P3_H
#define LARCONDTPCNV_LARcaliWAVEsubsetCNV_P3_H

#include "LArCondTPCnv/LArCaliWaveSubset_p3.h"
#include "LArRawConditions/LArCaliWave.h"
#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"
#include "LArRawConditions/LArConditionsSubset.h"

class MsgStream;

typedef LArConditionsSubset<LArCaliWaveVec>   LArCWTransType;  

class LArCaliWaveSubsetCnv_p3  
    : public T_AthenaPoolTPCnvConstBase<LArCWTransType, LArCaliWaveSubset_p3 >
{
public:
    using base_class::transToPers;
    using base_class::persToTrans;

    LArCaliWaveSubsetCnv_p3() {}
    virtual void   persToTrans(const LArCaliWaveSubset_p3*   persObj, LArCWTransType* transObj, MsgStream &log) const override;
    virtual void   transToPers(const LArCWTransType* transObj, LArCaliWaveSubset_p3*   persObj, MsgStream &log) const override;

};

#endif
