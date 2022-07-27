/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGMUONMONITORING_EFMUONMON_H
#define TRIGMUONMONITORING_EFMUONMON_H

#include "TrigMuonMonitorAlgorithm.h"
#include "xAODMuon/MuonContainer.h"

/*
This is a class for monitoring EFMuon.
 */
class EFMuonMon : public TrigMuonMonitorAlgorithm{

 public:
  EFMuonMon(const std::string& name, ISvcLocator* pSvcLocator );

  virtual StatusCode initialize() override;

 protected:
  virtual StatusCode fillVariablesPerChain(const EventContext &ctx, const std::string &chain) const override;
  virtual StatusCode fillVariablesPerOfflineMuonPerChain(const EventContext &ctx, const xAOD::Muon* mu, const std::string &chain) const override;
  virtual StatusCode fillVariables(const EventContext& ctx) const override;
  virtual StatusCode fillVariablesPerOfflineMuon(const EventContext &ctx, const xAOD::Muon* mu) const override;

 private:
  SG::ReadHandleKey<xAOD::MuonContainer> m_EFSAMuonContainerKey {this, "EFSAMuonContainerName", "HLT_Muons_RoI", "EFSAMuon container"};
  SG::ReadHandleKey<xAOD::MuonContainer> m_EFCBMuonContainerKey {this, "EFCBMuonContainerName", "HLT_MuonsCB_RoI", "EFCBMuon container"};
  SG::ReadHandleKey<xAOD::TrackParticleContainer> m_MStrackContainerKey {this, "ExtrapolatedMStrackConntainner", "HLT_MSExtrapolatedMuons_RoITrackParticles", "EFCBMuon container"};
  SG::ReadHandleKey<xAOD::TrackParticleContainer> m_CBtrackContainerKey {this, "CBtrackContainerName", "HLT_CBCombinedMuon_RoITrackParticles", "EFCBMuon container"};
  SG::ReadDecorHandleKey<xAOD::MuonContainer> m_muonIso30Key {this, "MuonIso03Name", "HLT_MuonsIso.ptcone03", "Isolation in ptcone03" };

  std::map<std::string, bool> m_doEFSA {};
  std::map<std::string, bool> m_doEFCB {};
  std::map<std::string, bool> m_doEFIso {};
};

#endif //TRIGMUONMONITORING_EFMUONMON_H
