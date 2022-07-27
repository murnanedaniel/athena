/**
 * Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration.
 *
 * @file HGTD_RawData/src/HGTD_RDO.cxx
 * @author Alexander Leopold <alexander.leopold@cern.ch>
 *
 * @brief Implementation of HGTD_RDO.h
 */

#include "HGTD_RawData/HGTD_RDO.h"

HGTD_RDO::HGTD_RDO(const Identifier rdo_id, const float toa,
                   const unsigned int tot, const unsigned short bcid,
                   const unsigned short l1_id,
                   const unsigned short l1_a)
    : Identifiable(),
      m_rdo_id(rdo_id),
      m_toa(toa),
      m_tot(tot),
      m_bcid(bcid),
      m_l1_id(l1_id),
      m_l1_a(l1_a) {}
