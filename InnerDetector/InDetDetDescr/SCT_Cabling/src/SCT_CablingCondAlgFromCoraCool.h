// -*- C++ -*-

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SCT_CablingCondAlgFromCoraCool_H
#define SCT_CablingCondAlgFromCoraCool_H
/**   
 *   @file SCT_CablingCondAlgFromCoraCool.h
 *
 *   @brief Fills an SCT cabling object from the database
 *
 *   @author Susumu Oda
 *   @date 08/03/2018
 */

//Athena includes
#include "AthenaBaseComps/AthReentrantAlgorithm.h"
#include "AthenaPoolUtilities/CondAttrListVec.h"
#include "SCT_Cabling/SCT_CablingData.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"

//Gaudi includes
#include "GaudiKernel/ServiceHandle.h"

//STL includes
#include <string>

//Forward declarations
class SCT_ID;

/**
 *    @class SCT_CablingCondAlgFromCoraCool
 *    @brief Condition algorithm which fills the SCT Cabling from the database, using CoraCool.
 *
 */

class SCT_CablingCondAlgFromCoraCool: public AthReentrantAlgorithm {
 public:

  SCT_CablingCondAlgFromCoraCool(const std::string& name, ISvcLocator* svc);
  virtual ~SCT_CablingCondAlgFromCoraCool() = default;
  virtual StatusCode initialize() override;
  virtual StatusCode execute(const EventContext& ctx) const override;
  virtual StatusCode finalize() override;
  /** Make this algorithm clonable. */
  virtual bool isClonable() const override { return true; };
  
private:

  bool insert(const IdentifierHash& hash, const SCT_OnlineId& onlineId, const SCT_SerialNumber& sn, SCT_CablingData* data) const;

  /// Input conditions folders
  SG::ReadCondHandleKey<CondAttrListVec> m_readKeyRod{this, "ReadKeyRod", "/SCT/DAQ/Config/ROD", "Key of input (raw) conditions folder of Rods"};
  SG::ReadCondHandleKey<CondAttrListVec> m_readKeyRodMur{this, "ReadKeyRodMur", "/SCT/DAQ/Config/RODMUR", "Key of input (raw) conditions folder of RodMurs"};
  SG::ReadCondHandleKey<CondAttrListVec> m_readKeyMur{this, "ReadKeyMur", "/SCT/DAQ/Config/MUR", "Key of input (raw) conditions folder of Murs"};
  SG::ReadCondHandleKey<CondAttrListVec> m_readKeyGeo{this, "ReadKeyGeo", "/SCT/DAQ/Config/Geog", "Key of input (raw) conditions folder of Geography"};
  /// Output condition data
  SG::WriteCondHandleKey<SCT_CablingData> m_writeKey{this, "WriteKey", "SCT_CablingData", "Key of output (derived) conditions data"};

  const SCT_ID* m_idHelper{nullptr};
};//end of class

#endif // SCT_CablingCondAlgFromCoraCool_H
