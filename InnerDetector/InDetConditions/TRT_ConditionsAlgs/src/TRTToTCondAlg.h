/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRTTOTCONDALG_H
#define TRTTOTCONDALG_H

#include <map>
#include <string>
#include <vector>

#include "AthenaBaseComps/AthAlgorithm.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "AthenaPoolUtilities/CondAttrListVec.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "Gaudi/Property.h"
#include "TRT_ConditionsData/TRTDedxcorrection.h"

class TRTToTCondAlg : public AthAlgorithm
{
 public:
  TRTToTCondAlg(const std::string& name, ISvcLocator* pSvcLocator);
  virtual ~TRTToTCondAlg() override;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;
  virtual StatusCode finalize() override;
  enum EDataBaseType {kOldDB,kNewDB,kNewDBOccCorr};
  StatusCode update1(TRTDedxcorrection& Dedxcorrection, const CondAttrListVec* channel_values);
  StatusCode update2(TRTDedxcorrection& Dedxcorrection, const CondAttrListCollection* attrListColl );
 
 protected:
  static void updateOldDBParameters(TRTDedxcorrection& Dedxcollection, std::map<std::string,std::vector<float> > &result_dict) ;
  static void updateNewDBParameters(TRTDedxcorrection& Dedxcorrection, std::map<std::string,std::vector<float> > &result_dict) ;
  static void updateOccupancyCorrectionParameters(TRTDedxcorrection & Dedxcorrection, std::map<std::string,std::vector<float> > &result_dict) ;

 private:
  SG::ReadCondHandleKey<CondAttrListVec> m_VecReadKey{this,"ToTVecReadKey","/TRT/Calib/ToT/ToTVectors","ToTVec in-key"};
  SG::ReadCondHandleKey<CondAttrListCollection> m_ValReadKey{this,"ToTValReadKey","/TRT/Calib/ToT/ToTValue","ToTVal in-key"};
  SG::WriteCondHandleKey<TRTDedxcorrection> m_WriteKey{this,"ToTWriteKey","Dedxcorrection","Dedxcorrection out-key"};

  static const std::vector<std::string> m_dictNamesOldDB;
  static const std::vector<std::string> m_dictNamesNewDB;
};
#endif
