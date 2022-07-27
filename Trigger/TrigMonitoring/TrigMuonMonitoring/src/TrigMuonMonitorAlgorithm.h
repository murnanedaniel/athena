/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGMUONMONITORING_TRIGMUONMONITORALGORITHM_H
#define TRIGMUONMONITORING_TRIGMUONMONITORALGORITHM_H

#include "AthenaMonitoring/AthMonitorAlgorithm.h"
#include "xAODMuon/MuonContainer.h"
#include "StoreGate/ReadHandleKey.h"
#include "MuonMatchingTool.h"


/**
 * @brief Base class from which analyzers can define a derived class to do specific analysis.
 * e.g. L2MuonSAMon is such a class already defined.
 * Analyzers should define @c fillVariablesXXX functions and these functions are automatically called in the @c fillHistograms function.
 * Analyzers can change @c selectEvents and @c selectMuons in their specific classes by overriding them if needed.
 * @todo Support monitoring algorithms using truth muons
 */
class TrigMuonMonitorAlgorithm : public AthMonitorAlgorithm {

 public:
  TrigMuonMonitorAlgorithm(const std::string& name, ISvcLocator* pSvcLocator );


  virtual StatusCode initialize() override;

  /**
   * @brief Function that steers anlayses.
   * It currently calles four types of analyses, @c fillVariables, @c fillVariablesPerOfflineMuon, @c fillVariablesPerChain and @c fillVariablesPerOfflineMuonPerChain
   * that can be overridden in subclasses to do specific analyses.
   * @see @c fillVariables, @c fillVariablesPerOfflineMuon, @c fillVariablesPerChain and @c fillVariablesPerOfflineMuonPerChain
   * @param ctx @c EventContext provided by athenaMT
   */
  virtual StatusCode fillHistograms(const EventContext &ctx) const override;

 protected:

  /**
   * @brief Function that defines the event selection for anlayses
   * User should reimlement in a subclass if needed.
   * @return True if the event is used for an analysis.
   */
  virtual bool selectEvents() const;

  /**
   * @brief Function that defines the event selection for anlayses
   * Users should reimlement in a subclass if needed.
   * @param muons Offline muons in the MuonContainer
   * @param probes List of offline muons that are used in analyses
   */
  virtual StatusCode selectMuons(SG::ReadHandle<xAOD::MuonContainer> &muons, std::vector<const xAOD::Muon*> &probes) const;

  /**
   * @brief Function that fills variables by just retrieving containers of trigger objects.
   * Users should reimlement in a subclass if needed.
   * @see @c fillHistograms
   * @param ctx @c EventContext provided by athenaMT
   */
  virtual StatusCode fillVariables(const EventContext &ctx) const;

  /**
   * @brief Function that fills variables that are compared to offline muons but the trigger chains are not specified.
   * This is called in the for loop of offline muons in @c fillHistograms.
   * Users should reimlement in a subclass if needed.
   * @see @c fillHistograms
   * @param ctx @c EventContext provided by athenaMT
   * @param mu Pointer to an offline muon provided in @c fillHistograms
   */
  virtual StatusCode fillVariablesPerOfflineMuon(const EventContext &ctx, const xAOD::Muon* mu) const;

  /**
   * @brief Function that fills variables of trigger objects associated to specified trigger chains.
   * This is called in the for loop of trigger chains in @c fillHistograms.
   * Users should reimlement in a subclass if needed.
   * @see @c fillHistograms
   * @param ctx @c EventContext provided by athenaMT
   * @param chain Trigger chain provided in @cfillHistograms
   */
  virtual StatusCode fillVariablesPerChain(const EventContext &ctx, const std::string &chain) const;

  /**
   * @brief Function that fills variables of trigger objects associated to specified trigger chains comparing offline muons.
   * This is called in the for loop of trigger chains and offline muons in @c fillHistograms.
   * Users should reimlement in a subclass if needed.
   * @see @c fillHistograms
   * @param ctx @c EventContext provided by athenaMT
   * @param mu Pointer to an offline muon provided in @c fillHistograms
   * @param chain Trigger chain provided in @c fillHistograms
   */
  virtual StatusCode fillVariablesPerOfflineMuonPerChain(const EventContext &ctx, const xAOD::Muon* mu, const std::string &chain) const;


  /**
   * @brief Function that fills variables of ratio plots.
   * @see @c TrigMuonMonitorAlgorithm.icc for the implementation
   * @param ctx @c EventContext provided by athenaMT
   * @param mu Pointer to an offline muon provided in @c fillHistograms
   * @param trigstep trigger step
   * @param type xAOD::Muon::TrackParticleType of offline muon
   * @param matchFunc Function pointer that implements cuts for the online muon candidates gotten by ReadHandle. 
   */
  template <class T, class FUNCT>
  StatusCode fillVariablesRatioPlots(const EventContext &ctx, const xAOD::Muon* mu,
                                     std::string &&trigstep,
                                     xAOD::Muon::TrackParticleType type,
                                     FUNCT matchFunc) const;

  /**
   * @brief Function that fills variables of etaphi2D plots.
   * @see @c TrigMuonMonitorAlgorithm.icc for the implementation
   * @param ctx @c EventContext provided by athenaMT
   * @param ReadHandleKey SG::ReadHandleKey of online muon.
   * @param trigstep trigger step
   * @param PosFunc Function pointer that implements cuts for the online muon candidates. 
   */
  template<class T>
  StatusCode fillVariableEtaPhi(const EventContext &ctx,
                                SG::ReadHandleKey<DataVector<T> > ReadHandleKey,
                                std::string &&trigstep,
                                std::tuple<bool,double,double> (*PosFunc)(const T*) = &TrigMuonMonitorAlgorithm::defaultPosFunc<T>) const;

  template<class T> static inline std::tuple<bool, double, double> defaultPosFunc(const T* trig);


  // ToolHandle
  ToolHandle<MuonMatchingTool> m_matchTool {this, "MuonMatchingTool", "MuonMatchingTool", "Tool for matching offline and online objects"};

  // ReadHandles
  SG::ReadHandleKey<xAOD::MuonContainer> m_MuonContainerKey {this, "MuonContainerName", "Muons", "Offline muon container"};

  // Properties
  /// List of trigger chains that are monitored in @c fillVariablesPerChain and @c fillVariablesPerOfflineMuonPerChain
  Gaudi::Property<std::vector<std::string> > m_monitored_chains {this, "MonitoredChains", {}, "Trigger chains that are monitored"};
  /// Requirement for the offline muon type considered in analyses
  Gaudi::Property<int> m_muontype {this, "MuonType", xAOD::Muon::MuonType::Combined, "MuonType used for monitoring"};
  /// Name of monitored group
  Gaudi::Property<std::string> m_group {this, "Group", "", "Histogram group"};

  /// Threshold for ratio measurement
  const float m_ratio_measurement_threshold = 4;

};

#include "TrigMuonMonitorAlgorithm.icc"

#endif //TRIGMUONMONITORING_TRIGMUONMONITORALGORITHM_H
