/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EVENTTAGTPCNV_DICT_H
#define EVENTTAGTPCNV_DICT_H
/**
 * @file EventTagTPCnvDict.h
 *
 * @brief Header file for dictionary generation
 *
 * @author  <Marcin.Nowak@cern.ch>
 */

#include "EventTagTPCnv/RawInfoSummaryForTag_p1.h"
#include "EventTagTPCnv/RawInfoSummaryForTagCnv_p1.h"

struct dummy {
  // template instances go here
  T_TPCnv<RawInfoSummaryForTag, RawInfoSummaryForTag_p1> m_conv1;
};

#endif

