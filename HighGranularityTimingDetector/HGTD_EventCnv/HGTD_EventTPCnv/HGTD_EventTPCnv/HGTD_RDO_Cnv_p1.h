/**
 * Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration.
 *
 * @file HGTD_EventTPCnv/HGTD_RDO_Cnv_p1.h
 * @author Alexander Leopold <alexander.leopold@cern.ch>
 * @brief Transient/Persistent converter for the HGTD_RDO class.
 *
 */

#ifndef HGTD_EVENTTPCNV_HGTD_RDO_CNV_P1_H
#define HGTD_EVENTTPCNV_HGTD_RDO_CNV_P1_H

#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"
#include "HGTD_EventTPCnv/HGTD_RDO_p1.h"
#include "HGTD_RawData/HGTD_RDO.h"

class MsgStream;

class HGTD_RDO_Cnv_p1
    : public T_AthenaPoolTPCnvBase<HGTD_RDO, HGTD_RDO_p1> {
public:
  HGTD_RDO_Cnv_p1() = default;

  void persToTrans(const HGTD_RDO_p1* pers_obj, HGTD_RDO* trans_obj,
                   MsgStream& log);

  void transToPers(const HGTD_RDO* trans_obj, HGTD_RDO_p1* pers_obj,
                   MsgStream& log);
};

#endif // HGTD_EVENTTPCNV_HGTD_RDO_CNV_P1_H
