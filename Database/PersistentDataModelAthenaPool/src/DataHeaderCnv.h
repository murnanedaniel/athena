/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DATAHEADERCNV_H
#define DATAHEADERCNV_H

/**
 *  @file DataHeaderCnv.h
 *  @brief This file contains the class definition for the DataHeaderCnv class.
 *  @author Peter van Gemmeren <gemmeren@anl.gov>
 **/

#include "AthenaPoolCnvSvc/T_AthenaPoolCustomCnv.h"
#include "PersistentDataModel/DataHeader.h"
#include "PersistentDataModelTPCnv/DataHeader_p6.h"
#include "PersistentDataModelTPCnv/DataHeaderCnv_p5.h"
#include "PersistentDataModelTPCnv/DataHeaderCnv_p6.h"

#include "GaudiKernel/IIncidentListener.h"

#include <map>

class  DataHeader_p5;
class  DataHeaderForm_p6;

/** @class DataHeaderCnv
 *  @brief This class provides the converter to customize the saving of DataHeader.
 **/
typedef  T_AthenaPoolCustomCnv<DataHeader, DataHeader_p6> DataHeaderCnvBase;
class DataHeaderCnv : public DataHeaderCnvBase, virtual public IIncidentListener
{
   friend class CnvFactory<DataHeaderCnv>;
public:
   DataHeaderCnv(ISvcLocator* svcloc);
   ~DataHeaderCnv();
   virtual StatusCode initialize() override;

   /// Extend base-class conversion methods
   virtual StatusCode updateRep(IOpaqueAddress* pAddress, DataObject* pObject) override;

   virtual StatusCode updateRepRefs(IOpaqueAddress* pAddress, DataObject* pObject) override;

   virtual StatusCode DataObjectToPool(DataObject* pObj, const std::string& tname) override;

   std::unique_ptr<DataHeader_p5> poolReadObject_p5();
   std::unique_ptr<DataHeader_p6> poolReadObject_p6();

   DataHeader_p6* createPersistent(DataHeader* transObj, DataHeaderForm_p6*) ;
   virtual DataHeader* createTransient() override;

   [[deprecated("this converter uses createPersistent() with 2 arguments")]]
   virtual DataHeader_p6* createPersistent(DataHeader*) override { return nullptr; }

   /// Incident service handle listening for EndInputFile.
   virtual void handle(const Incident& incident) override;
 
protected:
   DataHeaderCnv_p6     m_tpOutConverter;
   DataHeaderCnv_p6     m_tpInConverter;
  
   RootType             m_dhFormType;

   /// DHForm cache indexed by filename or reference for writing
  std::map<std::string,  std::unique_ptr<DataHeaderForm_p6> >    m_persFormMap;

   /// for use when reading DataHeader_p5
   DataHeaderCnv_p5                     m_tpInConverter_p5;  
   std::unique_ptr<DataHeaderForm_p5>   m_dhInForm5;
   std::string                          m_dhFormMdx;

   /// cached values for use with SharedWriter server
   DataHeader_p6*                       m_sharedWriterCachedDH = nullptr;  // no ownership
   /// map of cached DHForms for DataHeader ID
   std::map< std::string, std::unique_ptr<DataHeaderForm_p6> >  m_sharedWriterCachedDHForm;
   std::string                          m_sharedWriterCachedDHKey;
   std::string                          m_sharedWriterCachedDHToken;
};

#endif
