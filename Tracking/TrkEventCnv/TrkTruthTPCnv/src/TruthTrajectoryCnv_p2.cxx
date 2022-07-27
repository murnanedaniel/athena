/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// T/P converter for TruthTrajectory.
// Andrei Gaponenko <agaponenko@lbl.gov>, 2007
// Olivier Arnaez <olivier.arnaez@cern.ch>, 2015

#include "TrkTruthTPCnv/TruthTrajectoryCnv_p2.h"
#include "TrkTruthTPCnv/TruthTrajectory_p2.h"
#include "TrkTruthData/TruthTrajectory.h"

#include "GeneratorObjectsTPCnv/HepMcParticleLinkCnv_p2.h"

namespace {
  const HepMcParticleLinkCnv_p2 particleLinkConverter;
}

void TruthTrajectoryCnv_p2::persToTrans( const Trk::TruthTrajectory_p2* pers,
                                         TruthTrajectory* trans,
                                         MsgStream& msg ) const
{
  trans->resize(pers->size());
  for(Trk::TruthTrajectory_p2::size_type i=0; i<trans->size(); i++) {
    particleLinkConverter.persToTrans( &((*pers)[i]), &((*trans)[i]), msg);
  }
}

void TruthTrajectoryCnv_p2::transToPers( const TruthTrajectory* trans,
                                         Trk::TruthTrajectory_p2* pers,
                                         MsgStream& msg ) const
{
  pers->resize(trans->size());
  for(TruthTrajectory::size_type i=0; i<trans->size(); i++) {
    particleLinkConverter.transToPers( &((*trans)[i]), &((*pers)[i]), msg);
  }
}
