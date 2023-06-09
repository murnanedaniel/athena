/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DATAMODELATHENAPOOL_ELEMENTLINKCNV_P3_H
#define DATAMODELATHENAPOOL_ELEMENTLINKCNV_P3_H

/**
 *  @file ElementLinkCnv_p3.h
 *  @brief This file contains the class definition for the ElementLinkCnv_p3 class.
 *  @author Marcin.Nowak@cern.ch
 **/

#include "AthenaPoolCnvSvc/T_AthenaPoolTPConverter.h"
#include "ElementLink_p3.h"

/** @class ElementLinkCnv_p3<LINK>
 *  @brief Converter template for converters between transient ElementLink and its persistent representation. Template argument LINK is the type of the transient ElementLink. The type of the persistent link representation is automatically deducted - it can be integer index type, or string index type. More types can be added in ElementLink_p3.h
 **/

template <class LINK_TYPE>
class ElementLinkCnv_p3
   : public T_AthenaPoolTPCnvBase< LINK_TYPE, typename GeneratePersELinkType_p3<LINK_TYPE >::type > {
public:
  typedef	LINK_TYPE				Link_t;
  typedef 	typename GeneratePersELinkType_p3< Link_t >::type	PersLink_t;

  ElementLinkCnv_p3() {}
  
  virtual void persToTrans(const PersLink_t* pers, Link_t* trans, MsgStream& log) ;
  virtual void transToPers(const Link_t* trans, PersLink_t* pers, MsgStream& log) ;

  virtual void persToTrans(const PersLink_t& pers, Link_t& trans, MsgStream& log) ;
  virtual void transToPers(const Link_t& trans, PersLink_t& pers, MsgStream& log) ;

};

#ifndef __REFLEX__
#include "ElementLinkCnv_p3.icc"
#endif // not __REFLEX__

#endif
