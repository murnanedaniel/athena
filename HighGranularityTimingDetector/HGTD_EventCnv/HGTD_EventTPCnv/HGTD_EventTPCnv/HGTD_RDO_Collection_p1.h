/**
 * Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration.
 *
 * @file HGTD_EventTPCnv/HGTD_RDO_Collection_p1.h
 * @author Alexander Leopold <alexander.leopold@cern.ch>
 * @brief
 */

#ifndef HGTD_EVENTTPCNV_HGTD_RDO_COLLECTION_P1_H
#define HGTD_EVENTTPCNV_HGTD_RDO_COLLECTION_P1_H

#include "Identifier/IdentifierHash.h"

class HGTD_RDO_Collection_p1 {

public:
  typedef IdentifierHash::value_type IdType_t;

  HGTD_RDO_Collection_p1() : m_hash_id(0), m_size(0){};

  // container cnvs do conversions
  friend class HGTD_RDO_ContainerCnv_p1;

private:
  // Hash Identifier of this collection
  IdType_t m_hash_id;

  //  Number of RDOs in a given collection
  unsigned short m_size;
};

#endif // HGTD_EVENTTPCNV_HGTD_RDO_COLLECTION_P1_H
