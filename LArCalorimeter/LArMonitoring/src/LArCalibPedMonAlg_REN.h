/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef LARCALIBPEDMONALGREN_H
#define LARCALIBPEDMONALGREN_H

#include "AthenaMonitoring/AthMonitorAlgorithm.h"
#include "AthenaMonitoringKernel/Monitored.h"

#include "StoreGate/ReadHandleKey.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "LArRawEvent/LArCalibDigitContainer.h"
#include "LArRawEvent/LArAccumulatedDigitContainer.h"
#include "LArRawEvent/LArAccumulatedCalibDigitContainer.h"

#include "LArRawEvent/LArFebHeaderContainer.h"
#include "xAODEventInfo/EventInfo.h"

#include <string>
#include <vector>

class LArCalibPedMonAlgREN: public AthMonitorAlgorithm
{
 public:
  LArCalibPedMonAlgREN(const std::string& name,ISvcLocator* pSvcLocator );		      

  /** @brief Default destructor */
  virtual ~LArCalibPedMonAlgREN();

  /** @brief Overwrite dummy method from AlgTool */
  virtual StatusCode initialize() override;


  /** Called each event */
  virtual StatusCode fillHistograms( const EventContext& ctx ) const override;

 private:

  // keys to access info
  SG::ReadHandleKey<LArCalibDigitContainer> m_calibDigitContainerKey{this,"LArCalibDigitContainerKey","","SG key of LArCalibDigitContainer read from Bytestream"};
  SG::ReadHandleKey<LArAccumulatedDigitContainer> m_accDigitContainerKey{this,"LArAccumulatedDigitContainerKey","","SG key of LArAccumulatedDigitContainer read from Bytestream"};
  SG::ReadHandleKey<LArAccumulatedCalibDigitContainer> m_accCalibDigitContainerKey{this,"LArAccumulatedCalibDigitContainerKey","","SG key of LArAccumulatedCalibDigitContainer read from Bytestream"};
  SG::ReadHandleKey<LArFebHeaderContainer> m_hdrContKey{this, "LArFebHeaderKey", "LArFebHeader"};
  
  // Properties
  Gaudi::Property<std::vector<std::string> > m_partitions {this, "PartitionNames", {} };
  Gaudi::Property<std::vector<std::string> > m_SubDetNames{this, "SubDetNames", {} };
    //MonGroup(s) name
  Gaudi::Property<std::string> m_MonGroupName {this,"LArPedGroupName","LArPedMonGroup"};
  
    /* Histogram grouping (part) */
  std::vector<std::map<std::string,int> > m_histoGroups;

  mutable std::atomic<int>  m_nbOfFebBlocksTotal;
};

#endif

