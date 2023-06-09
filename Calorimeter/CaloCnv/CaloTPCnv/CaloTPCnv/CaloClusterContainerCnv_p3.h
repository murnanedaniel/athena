//Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef CALOTPCNV_CALOCLUSTERCONTAINERCNV_P3_H
#define CALOTPCNV_CALOCLUSTERCONTAINERCNV_P3_H

#include "CaloEvent/CaloClusterContainer.h"
#include "CaloEvent/CaloSamplingData.h"
#include "AthenaKernel/ITPCnvBase.h"
#include "CaloTPCnv/CaloTowerSegCnv_p1.h"
//#include "CaloTPCnv/CaloClusterMomentContainerCnv_p2.h"
#include "CaloTPCnv/CaloSamplingDataContainerCnv_p1.h"
#include "EventCommonTPCnv/P4EEtaPhiMCnv_p1.h"
#include "CaloEvent/CaloShowerContainer.h"
#include "CaloEvent/CaloCellLinkContainer.h"
#include "CaloTPCnv/CaloClusterContainer_p3.h"

#include "DataModelAthenaPool/ElementLinkCnv_p2.h"
#include "DataModel/ElementLink.h"
#include "AthenaPoolCnvSvc/ITPConverter.h"

class CaloClusterContainer;
class CaloCluster;

class CaloClusterContainerCnv_p3 : public ITPCnvBase {
public:
  CaloClusterContainerCnv_p3() {};
  virtual ~CaloClusterContainerCnv_p3() {}; 

  // Methods for invoking conversions on objects given by generic pointers.
  virtual void persToTransUntyped(const void* pers,
                                  void* trans,
                                  MsgStream& log);
  virtual void transToPersUntyped(const void* trans,
                                  void* pers,
                                  MsgStream& log);
  virtual const std::type_info& transientTInfo() const;

  /** return C++ type id of the persistent class this converter is for
      @return std::type_info&
  */
  virtual const std::type_info& persistentTInfo() const;

  void persToTrans(const CaloClusterContainer_p3*, CaloClusterContainer*, MsgStream &log);
  void transToPers(const CaloClusterContainer*, CaloClusterContainer_p3*, MsgStream &log);

private:
  //Conversion function for individual clusters (called in a loop over the container)
  void persToTrans(const CaloClusterContainer_p3::CaloCluster_p*, CaloCluster*, MsgStream &);
  void transToPers(const CaloCluster*, CaloClusterContainer_p3::CaloCluster_p*, MsgStream &);

  //Sub-Converters:
  CaloTowerSegCnv_p1                                     m_caloTowerSegCnv;
  P4EEtaPhiMCnv_p1                                       m_P4EEtaPhiMCnv;
//  CaloClusterMomentContainerCnv_p2                       m_momentContainerCnv;
  CaloSamplingDataContainerCnv_p1                        m_samplingDataContainerCnv;
  ElementLinkCnv_p2<ElementLink<CaloShowerContainer> >   m_showerElementLinkCnv;
  ElementLinkCnv_p2<ElementLink<CaloCellLinkContainer> > m_cellElementLinkCnv;
};


#endif
