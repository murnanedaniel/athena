/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonRDO/CscRawData.h"
#include "CscRawDataCnv_p4.h"
#include "GaudiKernel/GaudiException.h"

void
CscRawDataCnv_p4::transToPers(const CscRawData* trans, CscRawData_p4* pers, MsgStream &) 
{
    pers->m_amps        = trans->samples();
    pers->m_address     = trans->address();
    pers->m_id          = trans->identify();
    pers->m_time        = trans->time();
    pers->m_rpuID       = trans->rpuID();
    pers->m_width       = trans->width();
    pers->m_isTimeComputed        = trans->isTimeComputed();
    pers->m_hashId        = trans->hashId();
}

void
CscRawDataCnv_p4::persToTrans(const CscRawData_p4* pers, CscRawData* trans, MsgStream &) 
{
   *trans = CscRawData (pers->m_amps,
                        pers->m_address,
                        pers->m_id,
                        pers->m_rpuID,
                        pers->m_width);
   trans->setHashID (pers->m_hashId);
   if (pers->m_isTimeComputed)
     trans->setTime (pers->m_time);
}
