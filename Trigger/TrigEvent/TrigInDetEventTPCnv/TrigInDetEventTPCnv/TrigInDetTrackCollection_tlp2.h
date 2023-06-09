/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**********************************************************************************
 * @Project: Trigger
 * @Package: TrigInDetEventTPCnv
 * @Class  : TrigInDetTrackCollection_tlp1
 *
 * @brief "top level" persistent partner for TrigInDetTrackCollection
 *
 * @author Andrew Hamilton  <Andrew.Hamilton@cern.ch>  - U. Geneva
 * @author Francesca Bucci  <f.bucci@cern.ch>          - U. Geneva
 *
 * File and Version Information:
 * $Id: TrigInDetTrackCollection_tlp1.h,v 1.2 2009-04-01 22:08:44 ilija@vukotic.me Exp $
 **********************************************************************************/
#ifndef TRIGINDETEVENTTPCNV_TRIGINDETTRACKCOLLECTION_TLP2_H
#define TRIGINDETEVENTTPCNV_TRIGINDETTRACKCOLLECTION_TLP2_H

#include "TrigInDetTrackCollection_p1.h"
#include "TrigInDetTrack_p3.h"
#include "TrigInDetTrackFitPar_p3.h"

class TrigInDetTrackCollection_tlp2
{
 public:
  
  TrigInDetTrackCollection_tlp2() {}
  friend class TrigInDetTrackCollectionCnv_tlp2;
  
 private:
  
  std::vector< TrigInDetTrackCollection_p1 >  m_trigInDetTrackCollections;
  std::vector< TrigInDetTrackFitPar_p3 >      m_trigInDetTrackFitPars;
  std::vector< TrigInDetTrack_p3 >            m_trigInDetTracks;

};

#endif
