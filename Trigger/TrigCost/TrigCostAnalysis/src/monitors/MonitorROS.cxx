/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "MonitorROS.h"
#include "../counters/CounterROS.h"

#include <algorithm>

MonitorROS::MonitorROS(const std::string& name, const MonitoredRange* parent)
  : MonitorBase(name, parent) {
}


StatusCode MonitorROS::newEvent(const CostData& data, const float weight) {

  // Prepare ROB id per corresponding ROS name map
  if (m_robToRos.empty()) {
    const std::map<std::string, std::vector<uint32_t>>& rosToRobMap = data.rosToRobMap();
    for (const auto& rosRequest : rosToRobMap) {
      for (uint32_t robId : rosRequest.second) {
        m_robToRos[robId] = rosRequest.first;
      }
    }
  }

  if (data.rosCollection().empty()){
    ATH_MSG_DEBUG("The ROS collection is empty!");
  }

  for (const xAOD::TrigComposite* tc : data.rosCollection()) {
    auto robIds = tc->getDetail<std::vector<uint32_t>>("robs_id");
    // Create set of unique ROS for this request
    std::set<std::string> rosPerRequest;
    for (uint32_t robId : robIds) {
      if (!m_robToRos.count(robId)){
        ATH_MSG_WARNING("ROS for ROB 0x" << std::hex << robId << " is missing");
        continue;
      }
      rosPerRequest.insert(m_robToRos.at(robId));
    }

    for (const std::string& rosName : rosPerRequest) {
      ATH_CHECK( getCounter(rosName)->newEvent(data, tc->index(), weight) );
    }
  }

  return StatusCode::SUCCESS;
}


std::unique_ptr<CounterBase> MonitorROS::newCounter(const std::string& name) {
  return std::make_unique<CounterROS>(name, this);
} 
