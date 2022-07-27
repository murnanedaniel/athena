/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/
#include "MuonEventTPCnv/MuonRDO/MM_RawDataCnv_p2.h"

void Muon::MM_RawDataCnv_p2::persToTrans( const Muon::MM_RawData_p2 *, Muon::MM_RawData *,MsgStream & /**log*/ ) 
{
  // Not much point calling this!
}

void Muon::MM_RawDataCnv_p2::transToPers( const Muon::MM_RawData *transObj, Muon::MM_RawData_p2 *persObj, MsgStream & /**log*/ )
{
  persObj->m_id = transObj->m_id.get_identifier32().get_compact();
  persObj->m_channel = transObj->channel();
  persObj->m_time = transObj->time();
  persObj->m_charge = transObj->charge();
  persObj->m_relBcid = transObj->relBcid();
}

Muon::MM_RawData* Muon::MM_RawDataCnv_p2::createTransient(const Muon::MM_RawData_p2* persObj, MsgStream& /**log*/)
{
// since the persistent class used uint for the member variables and the transient one uses int,
// the static cast is added here explicitly  
Muon::MM_RawData*  trans = new MM_RawData( Identifier (persObj->m_id),
					    static_cast<int>(persObj->m_channel),
					     static_cast<int>(persObj->m_time), 
					     static_cast<int>(persObj->m_charge),
               persObj->m_relBcid,
	       false	 ); // any RDOs produced with this persistent version have the time and the charge in physical units. pscholer March 2022 
  
  return trans;
}

