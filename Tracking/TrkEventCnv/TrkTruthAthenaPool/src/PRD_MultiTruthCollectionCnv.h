/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// Dear emacs, this is -*-c++-*-
#ifndef PRD_MULTITRUTHCOLLECTIONCNV_H
#define PRD_MULTITRUTHCOLLECTIONCNV_H

#include "TrkTruthData/PRD_MultiTruthCollection.h"
#include "TrkTruthTPCnv/PRD_MultiTruthCollectionCnv_p1.h"
#include "TrkTruthTPCnv/PRD_MultiTruthCollectionCnv_p2.h"
#include "TrkTruthTPCnv/PRD_MultiTruthCollectionCnv_p3.h"
#include "AthenaPoolCnvSvc/T_AthenaPoolCustomCnv.h"

namespace Trk { class PRD_MultiTruthCollection_p3; }

typedef Trk::PRD_MultiTruthCollection_p3 PRD_MultiTruthCollectionPERS;


typedef T_AthenaPoolCustomCnv<PRD_MultiTruthCollection,PRD_MultiTruthCollectionPERS> PRD_MultiTruthCollectionCnvBase;

class PRD_MultiTruthCollectionCnv : public PRD_MultiTruthCollectionCnvBase 
{
  friend class CnvFactory<PRD_MultiTruthCollectionCnv>;
protected:
public:
  PRD_MultiTruthCollectionCnv(ISvcLocator*);
protected:
  virtual PRD_MultiTruthCollection* createTransient();
  virtual PRD_MultiTruthCollectionPERS* createPersistent(PRD_MultiTruthCollection*);
private:
  PRD_MultiTruthCollectionCnv_p1 m_converter_p1;
  PRD_MultiTruthCollectionCnv_p2 m_converter_p2;
  PRD_MultiTruthCollectionCnv_p3 m_converter_p3;
  static const pool::Guid s_p0_guid;
  static const pool::Guid s_p1_guid;
  static const pool::Guid s_p2_guid;
  static const pool::Guid s_p3_guid;
};

#endif/*CALOATHENAPOOL_CALOCELLLINKCONTAINERCNV_H*/
