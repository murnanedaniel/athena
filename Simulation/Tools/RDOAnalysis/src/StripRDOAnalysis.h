/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/


#ifndef STRIP_RDO_ANALYSIS_H
#define STRIP_RDO_ANALYSIS_H

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ITHistSvc.h"
#include "StoreGate/ReadHandleKey.h"

#include "InDetRawData/InDetRawDataCLASS_DEF.h"
#include "InDetRawData/InDetRawDataContainer.h"
#include "InDetSimData/InDetSimDataCollection.h"
#include "InDetSimData/SCT_SimHelper.h"
#include "InDetIdentifier/SCT_ID.h"

#include <string>
#include <vector>
#include "TH1.h"
#include "TH2.h"

class TTree;
class SCT_ID;
class SCT_RDORawData;

namespace InDetDD {
   class SCT_DetectorManager;
}

namespace ITk
{

class StripRDOAnalysis : public AthAlgorithm
{

public:
  StripRDOAnalysis(const std::string& name, ISvcLocator* pSvcLocator);

  virtual StatusCode initialize() override final;
  virtual StatusCode execute() override final;

private:
  SG::ReadHandleKey<SCT_RDO_Container> m_inputKey {this, "CollectionName", "ITkStripRDOs", "Input ITk Strip RDO collection name"};
  SG::ReadHandleKey<InDetSimDataCollection> m_inputTruthKey {this, "SDOCollectionName", "ITkStripSDO_Map", "Input ITk Strip SDO collection name"};
  const SCT_ID *m_sctID{};
  const InDetDD::SCT_DetectorManager *m_SCT_Manager{};

  Gaudi::Property<std::string> m_histPath {this, "HistPath", "/RDOAnalysis/ITkStrip/", ""};
  Gaudi::Property<std::string> m_sharedHistPath {this, "SharedHistPath", "/RDOAnalysis/histos/", ""};
  Gaudi::Property<std::string> m_ntuplePath {this, "NtuplePath", "/RDOAnalysis/ntuples/", ""};
  Gaudi::Property<std::string> m_ntupleName {this, "NtupleName", "ITkStrip", ""};
  Gaudi::Property<bool> m_doPosition {this, "DoPosition", true, ""};

  ServiceHandle<ITHistSvc> m_thistSvc {this, "HistSvc", "THistSvc", ""};

  // RDO
  std::vector<unsigned long long>* m_rdoID{};
  std::vector<unsigned int>* m_rdoWord{};
  // SCT_ID
  std::vector<int>* m_barrelEndcap{};
  std::vector<int>* m_layerDisk{};
  std::vector<int>* m_phiModule{};
  std::vector<int>* m_etaModule{};
  std::vector<int>* m_side{};
  std::vector<int>* m_strip{};
  std::vector<int>* m_row{};
  // SCT_RDORawData
  std::vector<int>* m_groupSize{};
  // Global and Local positions
  std::vector<double>* m_globalX0{};
  std::vector<double>* m_globalY0{};
  std::vector<double>* m_globalZ0{};
  std::vector<double>* m_globalX1{};
  std::vector<double>* m_globalY1{};
  std::vector<double>* m_globalZ1{};
  std::vector<double>* m_localX{};
  std::vector<double>* m_localY{};
  std::vector<double>* m_localZ{};


  // SDO
  std::vector<unsigned long long>* m_sdoID{};
  std::vector<int>* m_sdoWord{};
  // SCT_ID
  std::vector<int>* m_barrelEndcap_sdo{};
  std::vector<int>* m_layerDisk_sdo{};
  std::vector<int>* m_phiModule_sdo{};
  std::vector<int>* m_etaModule_sdo{};
  std::vector<int>* m_side_sdo{};
  std::vector<int>* m_strip_sdo{};
  std::vector<int>* m_row_sdo{};
  // SCT_SimHelper
  std::vector<bool>* m_noise{};
  std::vector<bool>* m_belowThresh{};
  std::vector<bool>* m_disabled{};
  // Deposit - particle link + energy (charge)
  std::vector<int>* m_barcode{};
  std::vector<int>* m_eventIndex{};
  std::vector<float>* m_charge{};
  std::vector< std::vector<int> >* m_barcode_vec{};
  std::vector< std::vector<int> >* m_eventIndex_vec{};
  std::vector< std::vector<float> >* m_charge_vec{};

  // HISTOGRAMS
  TH1* m_h_rdoID{};
  TH1* m_h_rdoWord{};
  TH1* m_h_barrelEndcap{};
  TH1* m_h_layerDisk{};
  TH1* m_h_phiModule{};
  TH1* m_h_etaModule{};
  TH1* m_h_side{};
  TH1* m_h_strip{};
  TH1* m_h_row{};
  TH1* m_h_groupSize{};
  TH2* m_h_phi_v_eta{};
  // barrel SCT
  TH1* m_h_brlLayer{};
  TH1* m_h_brlPhiMod{};
  TH1* m_h_brlEtaMod{};
  TH1* m_h_brlSide{};
  TH1* m_h_brlStrip{};
  TH1* m_h_brlGroupSize{};
  TH2* m_h_brl_phi_v_eta{};
  // endcap SCT
  TH1* m_h_ecDisk{};
  TH1* m_h_ecPhiMod{};
  TH1* m_h_ecEtaMod{};
  TH1* m_h_ecSide{};
  TH1* m_h_ecStrip{};
  TH1* m_h_ecGroupSize{};
  TH2* m_h_ec_phi_v_eta{};

  TH1* m_h_sdoID{};
  TH1* m_h_sdoWord{};
  TH1* m_h_barrelEndcap_sdo{};
  TH1* m_h_layerDisk_sdo{};
  TH1* m_h_phiModule_sdo{};
  TH1* m_h_etaModule_sdo{};
  TH1* m_h_side_sdo{};
  TH1* m_h_strip_sdo{};
  TH1* m_h_row_sdo{};
  TH1* m_h_barcode{};
  TH1* m_h_eventIndex{};
  TH1* m_h_charge{};
  TH2* m_h_phi_v_eta_sdo{};

  TH1* m_h_belowThresh_brl{};
  TH1* m_h_belowThresh_ec{};

  TH1* m_h_disabled_brl{};
  TH1* m_h_disabled_ec{};

  std::vector<TH1*> m_h_brl_strip_perLayer;
  std::vector<TH1*> m_h_ec_strip_perLayer;

  TH2* m_h_globalXY{};
  TH2* m_h_globalZR{};
  LockedHandle<TH2> m_h_globalXY_shared{};
  LockedHandle<TH2> m_h_globalZR_shared{};
  TH1* m_h_globalX{};
  TH1* m_h_globalY{};
  TH1* m_h_globalZ{};

  TH1* m_h_truthMatchedRDOs{};

  TTree* m_tree{};

};

} // namespace ITk

#endif // STRIP_RDO_ANALYSIS_H
