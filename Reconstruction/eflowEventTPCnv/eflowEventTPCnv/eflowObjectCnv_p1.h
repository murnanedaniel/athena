/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFLOWOBJECTCNV_P1_h
#define EFLOWOBJECTCNV_P1_h


#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"
#include "DataModelAthenaPool/ElementLinkVectorCnv_p1.h"
#include "DataModelAthenaPool/ElementLinkCnv_p1.h"


#include "eflowEvent/eflowObject.h"
#include "eflowEventTPCnv/eflowObject_p1.h"

#include "VxVertex/VxContainer.h"
#include "CaloEvent/CaloClusterContainer.h"
#include "Particle/TrackParticleContainer.h"



class MsgStream;

class eflowObjectCnv_p1 : public T_AthenaPoolTPCnvBase<eflowObject, eflowObject_p1>
{

 public:

  eflowObjectCnv_p1 () {};

  virtual void persToTrans( const eflowObject_p1 *persObj, eflowObject *transObj, MsgStream &eflowStream);
  virtual void transToPers( const eflowObject *transObj,eflowObject_p1 *persObj, MsgStream &eflowStream );

  friend class eflowObject;
  friend class eflowObject_p1;

 private:

  ElementLinkCnv_p1<ElementLink<VxContainer> > m_elVxCnv;
  ElementLinkVectorCnv_p1<ElementLinkVector<Rec::TrackParticleContainer> > m_navTrackCnv;
  ElementLinkVectorCnv_p1<ElementLinkVector<CaloClusterContainer> > m_navClusCnv;

};
#endif
