/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "StripRDOAnalysis.h"
#include "SCT_ReadoutGeometry/SCT_DetectorManager.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "ReadoutGeometryBase/SiLocalPosition.h"
#include "GeneratorObjects/McEventCollection.h"
#include "StoreGate/ReadHandle.h"
#include "TTree.h"
#include "TString.h"

#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>

namespace ITk
{

StripRDOAnalysis::StripRDOAnalysis(const std::string& name, ISvcLocator *pSvcLocator)
  : AthAlgorithm(name, pSvcLocator)
{
}

StatusCode StripRDOAnalysis::initialize() {
  ATH_MSG_DEBUG( "Initializing SCT_RDOAnalysis" );

  // This will check that the properties were initialized
  // properly by job configuration.
  ATH_CHECK( m_inputKey.initialize() );
  ATH_CHECK( m_inputTruthKey.initialize() );

  // Grab SCT_ID helper
  ATH_CHECK(detStore()->retrieve(m_sctID, "SCT_ID"));

  ATH_CHECK(detStore()->retrieve(m_SCT_Manager, "ITkStrip"));
  // Grab Ntuple and histogramming service for tree
  ATH_CHECK(m_thistSvc.retrieve());

  m_tree = new TTree(m_ntupleName.value().c_str(), "ITkStripRDOAnalysis");
  ATH_CHECK(m_thistSvc->regTree(m_ntuplePath.value() + m_ntupleName.value(), m_tree));
  if (m_tree) {
    // ITk Strip RDO
    m_tree->Branch("rdoID", &m_rdoID);
    m_tree->Branch("rdoWord", &m_rdoWord);
    m_tree->Branch("barrelEndcap", &m_barrelEndcap);
    m_tree->Branch("layerDisk", &m_layerDisk);
    m_tree->Branch("phiModule", &m_phiModule);
    m_tree->Branch("etaModule", &m_etaModule);
    m_tree->Branch("side", &m_side);
    m_tree->Branch("strip", &m_strip);
    m_tree->Branch("row", &m_row);
    m_tree->Branch("groupSize", &m_groupSize);
    // Global coordinates
    if (m_doPosition) {
      m_tree->Branch("globalX0", &m_globalX0);
      m_tree->Branch("globalY0", &m_globalY0);
      m_tree->Branch("globalZ0", &m_globalZ0);
      m_tree->Branch("globalX1", &m_globalX1);
      m_tree->Branch("globalY1", &m_globalY1);
      m_tree->Branch("globalZ1", &m_globalZ1);
      m_tree->Branch("localX", &m_localX);
      m_tree->Branch("localY", &m_localY);
      m_tree->Branch("localZ", &m_localZ);
    }
    // ITk Strip SDO deposits
    m_tree->Branch("sdoID", &m_sdoID);
    m_tree->Branch("sdoWord", &m_sdoWord);
    m_tree->Branch("barrelEndcap_sdo", &m_barrelEndcap_sdo);
    m_tree->Branch("layerDisk_sdo", &m_layerDisk_sdo);
    m_tree->Branch("phiModule_sdo", &m_phiModule_sdo);
    m_tree->Branch("etaModule_sdo", &m_etaModule_sdo);
    m_tree->Branch("side_sdo", &m_side_sdo);
    m_tree->Branch("strip_sdo", &m_strip_sdo);
    m_tree->Branch("row_sdo", &m_row_sdo);
    m_tree->Branch("noise", &m_noise);
    m_tree->Branch("belowThresh", &m_belowThresh);
    m_tree->Branch("disabled", &m_disabled);
    m_tree->Branch("barcode", &m_barcode);
    m_tree->Branch("eventIndex", &m_eventIndex);
    m_tree->Branch("charge", &m_charge);
    m_tree->Branch("barcode_vec", &m_barcode_vec);
    m_tree->Branch("eventIndex_vec", &m_eventIndex_vec);
    m_tree->Branch("charge_vec", &m_charge_vec);
  }
  else {
    ATH_MSG_ERROR("No tree found!");
  }

  // HISTOGRAMS
  m_h_rdoID = new TH1F("h_rdoID", "rdoID", 100, 0, 25e17);
  m_h_rdoID->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_rdoID->GetName(), m_h_rdoID));


  m_h_rdoWord = new TH1F("h_rdoWord", "rdoWord", 100, 0, 17e6);
  m_h_rdoWord->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_rdoWord->GetName(), m_h_rdoWord));

  m_h_barrelEndcap = new TH1F("h_barrelEndcap", "Barrel or Endcap", 100, -3, 3);
  m_h_barrelEndcap->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_barrelEndcap->GetName(), m_h_barrelEndcap));

  m_h_layerDisk = new TH1F("h_layerDisk", "Barrel layer or Endcap disk", 100, 0, 10);
  m_h_layerDisk->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_layerDisk->GetName(), m_h_layerDisk));

  m_h_phiModule = new TH1F("h_phiModule", "Phi module", 100, 0, 60);
  m_h_phiModule->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_phiModule->GetName(), m_h_phiModule));

  m_h_etaModule = new TH1F("h_etaModule", "Eta module", 121, -60, 60);
  m_h_etaModule->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_etaModule->GetName(), m_h_etaModule));

  m_h_side = new TH1F("h_side", "Side", 100, 0, 1.5);
  m_h_side->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_side->GetName(), m_h_side));

  m_h_strip = new TH1F("h_strip", "Strip", 100, 0, 800);
  m_h_strip->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_strip->GetName(), m_h_strip));

  m_h_row = new TH1F("h_row", "Row", 100, 0, 4.5);
  m_h_row->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_row->GetName(), m_h_row));

  m_h_groupSize = new TH1F("h_groupSize", "Group size", 100, 0, 150);
  m_h_groupSize->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_groupSize->GetName(), m_h_groupSize));

  m_h_phi_v_eta = new TH2F("h_phi_v_eta", "Phi module vs eta module", 100, -7, 7, 100, 0, 60);
  m_h_phi_v_eta->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_phi_v_eta->GetName(), m_h_phi_v_eta));

  m_h_brlLayer = new TH1F("h_brlLayer", "Barrel layer", 100, 0, 10);
  m_h_brlLayer->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_brlLayer->GetName(), m_h_brlLayer));

  m_h_brlPhiMod = new TH1F("h_brlPhiMod", "Barrel phi module", 100, 0, 60);
  m_h_brlPhiMod->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_brlPhiMod->GetName(), m_h_brlPhiMod));

  m_h_brlEtaMod = new TH1F("h_brlEtaMod", "Barrel eta module", 121, -60, 60);
  m_h_brlEtaMod->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_brlEtaMod->GetName(), m_h_brlEtaMod));

  m_h_brlSide = new TH1F("h_brlSide", "Barrel side", 100, 0, 1.5);
  m_h_brlSide->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_brlSide->GetName(), m_h_brlSide));

  m_h_brlStrip = new TH1F("h_brlStrip", "Barrel strip", 100, 0, 800);
  m_h_brlStrip->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_brlStrip->GetName(), m_h_brlStrip));

  m_h_brlGroupSize = new TH1F("h_brlGroupSize", "Barrel group size", 100, 0, 150);
  m_h_brlGroupSize->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_brlGroupSize->GetName(), m_h_brlGroupSize));

  m_h_brl_phi_v_eta = new TH2F("h_brl_phi_v_eta", "Barrel phi module vs eta module", 100, -7, 7, 100, 0, 60);
  m_h_brl_phi_v_eta->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_brl_phi_v_eta->GetName(), m_h_brl_phi_v_eta));

  m_h_ecDisk = new TH1F("h_ecDisk", "Endcap disk", 100, 0, 10);
  m_h_ecDisk->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_ecDisk->GetName(), m_h_ecDisk));

  m_h_ecPhiMod = new TH1F("h_ecPhiMod", "Endcap phi module", 100, 0, 60);
  m_h_ecPhiMod->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_ecPhiMod->GetName(), m_h_ecPhiMod));

  m_h_ecEtaMod = new TH1F("h_ecEtaMod", "Endcap eta module", 21, 0, 20);
  m_h_ecEtaMod->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_ecEtaMod->GetName(), m_h_ecEtaMod));

  m_h_ecSide = new TH1F("h_ecSide", "Endcap side", 100, 0, 1.5);
  m_h_ecSide->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_ecSide->GetName(), m_h_ecSide));

  m_h_ecStrip = new TH1F("h_ecStrip", "Endcap strip", 100, 0, 800);
  m_h_ecStrip->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_ecStrip->GetName(), m_h_ecStrip));

  m_h_ecGroupSize = new TH1F("h_ecGroupSize", "Endcap group size", 100, 0, 150);
  m_h_ecGroupSize->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_ecGroupSize->GetName(), m_h_ecGroupSize));

  m_h_ec_phi_v_eta = new TH2F("h_ec_phi_v_eta", "Endcap phi module vs eta module", 100, -7.5, 7.5, 100, 0, 60);
  m_h_ec_phi_v_eta->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_ec_phi_v_eta->GetName(), m_h_ec_phi_v_eta));

  m_h_sdoID = new TH1F("h_sdoID", "sdoID", 100, 0, 1e18);
  m_h_sdoID->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_sdoID->GetName(), m_h_sdoID));

  m_h_sdoWord = new TH1F("h_sdoWord", "sdoWord", 100, 0, 1e7);
  m_h_sdoWord->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_sdoWord->GetName(), m_h_sdoWord));

  m_h_barrelEndcap_sdo = new TH1F("h_barrelEndcap_sdo", "Barrel or Endcap (SDO)", 100, -3, 3);
  m_h_barrelEndcap_sdo->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_barrelEndcap_sdo->GetName(), m_h_barrelEndcap_sdo));

  m_h_layerDisk_sdo = new TH1F("h_layerDisk_sdo", "Barrel layer or Endcap disk (SDO)", 100, 0, 10);
  m_h_layerDisk_sdo->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_layerDisk_sdo->GetName(), m_h_layerDisk_sdo));

  m_h_phiModule_sdo = new TH1F("h_phiModule_sdo", "Phi module (SDO)", 100, 0, 60);
  m_h_phiModule_sdo->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_phiModule_sdo->GetName(), m_h_phiModule_sdo));

  m_h_etaModule_sdo = new TH1F("h_etaModule_sdo", "Eta module (SDO)", 121, -60, 60);
  m_h_etaModule_sdo->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_etaModule_sdo->GetName(), m_h_etaModule_sdo));

  m_h_side_sdo = new TH1F("h_side_sdo", "Side (SDO)", 100, 0, 1.5);
  m_h_side_sdo->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_side_sdo->GetName(), m_h_side_sdo));

  m_h_strip_sdo = new TH1F("h_strip_sdo", "Strip (SDO)", 100, 0, 800);
  m_h_strip_sdo->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_strip_sdo->GetName(), m_h_strip_sdo));

  m_h_row_sdo = new TH1F("h_row_sdo", "Row (SDO)", 100, 0, 4.5);
  m_h_row_sdo->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_row_sdo->GetName(), m_h_row_sdo));

  m_h_barcode = new TH1F("h_barcode", "Barcode (SDO)", 100, 0, 2.2e5);
  m_h_barcode->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_barcode->GetName(), m_h_barcode));

  m_h_eventIndex = new TH1F("h_eventIndex", "Event index (SDO)", 100, 0, 10);
  m_h_eventIndex->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_eventIndex->GetName(), m_h_eventIndex));

  m_h_charge = new TH1F("h_charge", "Charge (SDO)", 100, 0, 6e6);
  m_h_charge->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_charge->GetName(), m_h_charge));

  m_h_phi_v_eta_sdo = new TH2F("h_phi_v_eta_sdo", "Phi module vs eta module (SDO)", 100, -7, 7, 100, 0, 60);
  m_h_phi_v_eta_sdo->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_phi_v_eta_sdo->GetName(), m_h_phi_v_eta_sdo));

  m_h_belowThresh_brl = new TH1F("h_belowThresh_brl", "Below threshold strips - Barrel; # below threshold strips; layer", 8, -0.5, 7.5);
  m_h_belowThresh_brl->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_belowThresh_brl->GetName(), m_h_belowThresh_brl));

  m_h_belowThresh_ec = new TH1F("h_belowThresh_ec", "Below threshold strips - Endcap; # below threshold strips; layer", 8, -0.5, 7.5);
  m_h_belowThresh_ec->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_belowThresh_ec->GetName(), m_h_belowThresh_ec));

  m_h_disabled_brl = new TH1F("m_h_disabled_brl", "Disabled strips - Barrel; # disabled strips; layer", 8, -0.5, 7.5);
  m_h_disabled_brl->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_disabled_brl->GetName(), m_h_disabled_brl));

  m_h_disabled_ec = new TH1F("m_h_disabled_ec", "Disabled strips - Endcap; # disabled strips; layer", 8, -0.5, 7.5);
  m_h_disabled_ec->StatOverflows();
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_disabled_ec->GetName(), m_h_disabled_ec));

  for (unsigned int layer=0; layer<4; layer++) {
    m_h_brl_strip_perLayer.emplace_back(new TH1F(("m_h_brl_strip_perLayer"+std::to_string(layer)).c_str(), ("Strip index - Barrel - Layer "+std::to_string(layer)).c_str(), 1300, 0, 1300));
    m_h_brl_strip_perLayer.back()->StatOverflows();
    ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_brl_strip_perLayer.back()->GetName(), m_h_brl_strip_perLayer.back()));
  }

  for (unsigned int layer=0; layer<9; layer++) {
    m_h_ec_strip_perLayer.emplace_back(new TH1F(("m_h_ec_strip_perLayer"+std::to_string(layer)).c_str(), ("Strip index - Barrel - Layer "+std::to_string(layer)).c_str(), 1300, 0, 1300));
    m_h_ec_strip_perLayer.back()->StatOverflows();
    ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_ec_strip_perLayer.back()->GetName(), m_h_ec_strip_perLayer.back()));
  }

  m_h_globalXY = new TH2F("m_h_globalXY","m_h_globalXY; x [mm]; y [mm]",2200,-1100.,1100.,2200,1100.,1100.);
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_globalXY->GetName(), m_h_globalXY));
  m_h_globalZR = new TH2F("m_h_globalZR","m_h_globalZR; z [mm]; r [mm]",6800,-3400.,3400.,1100,0.,1100.);
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_globalZR->GetName(), m_h_globalZR));
  m_h_globalX = new TH1F("m_h_globalX","m_h_globalX; x [mm]",2200,-1100.,1100.);
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_globalX->GetName(), m_h_globalX));
  m_h_globalY = new TH1F("m_h_globalY","m_h_globalY; y [mm]",2200,-1100.,1100.);
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_globalY->GetName(), m_h_globalY));
  m_h_globalZ = new TH1F("m_h_globalZ","m_h_globalZ; z [mm]",6800,-3400.,3400.);
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_globalZ->GetName(), m_h_globalZ));

  // Special shared ITk histograms
  std::string xy_name = "h_ITk_xy";
  auto xy = std::make_unique<TH2D>(xy_name.c_str(), xy_name.c_str(), 2200, -1100, 1100, 2200, -1100, 1100);
  xy->StatOverflows();
  ATH_CHECK(m_thistSvc->regShared(m_sharedHistPath + xy_name, std::move(xy), m_h_globalXY_shared));

  std::string zr_name = "h_ITk_zr";
  auto zr = std::make_unique<TH2D>(zr_name.c_str(), zr_name.c_str(), 6800, -3400, 3400, 1100, 0, 1100);
  zr->StatOverflows();
  ATH_CHECK(m_thistSvc->regShared(m_sharedHistPath + zr_name, std::move(zr), m_h_globalZR_shared));

  m_h_truthMatchedRDOs = new TH1F("h_TruthMatchedITkStripRDOs", "h_TruthMatchedITkStripRDOs", 4, 1, 5);
  TString truthMatchBinLables[4] = { "All RDOs", "Truth Matched", "HS Matched", "Unmatched" };
  for(unsigned int ibin = 1; ibin < 5; ibin++) {
    m_h_truthMatchedRDOs->GetXaxis()->SetBinLabel(ibin, truthMatchBinLables[ibin-1]);
  }
  ATH_CHECK(m_thistSvc->regHist(m_histPath + m_h_truthMatchedRDOs->GetName(), m_h_truthMatchedRDOs));

  return StatusCode::SUCCESS;
}

StatusCode StripRDOAnalysis::execute() {
  ATH_MSG_DEBUG( "In StripRDOAnalysis::execute()" );

  m_rdoID->clear();
  m_rdoWord->clear();
  m_barrelEndcap->clear();
  m_layerDisk->clear();
  m_phiModule->clear();
  m_etaModule->clear();
  m_side->clear();
  m_strip->clear();
  m_row->clear();
  m_groupSize->clear();
  if (m_doPosition) {
    m_globalX0->clear();
    m_globalY0->clear();
    m_globalZ0->clear();
    m_globalX1->clear();
    m_globalY1->clear();
    m_globalZ1->clear();
    m_localX->clear();
    m_localY->clear();
    m_localZ->clear();
  }
  m_sdoID->clear();
  m_sdoWord->clear();
  m_barrelEndcap_sdo->clear();
  m_layerDisk_sdo->clear();
  m_phiModule_sdo->clear();
  m_etaModule_sdo->clear();
  m_side_sdo->clear();
  m_strip_sdo->clear();
  m_row_sdo->clear();
  m_noise->clear();
  m_belowThresh->clear();
  m_disabled->clear();
  m_barcode->clear();
  m_eventIndex->clear();
  m_charge->clear();
  m_barcode_vec->clear();
  m_eventIndex_vec->clear();
  m_charge_vec->clear();

  // RawData
  SG::ReadHandle<SCT_RDO_Container> p_SCT_RDO_cont (m_inputKey);
  //Adding SimMap and McEvent here for added truthMatching checks
  SG::ReadHandle<InDetSimDataCollection> simDataMapSCT (m_inputTruthKey);
  SG::ReadHandle<McEventCollection> mcEventCollection("TruthEvent");

  const HepMC::GenEvent* hardScatterEvent(nullptr);
  bool doTruthMatching = true;
  if (mcEventCollection->size()==0){
    ATH_MSG_WARNING("Failed to retrieve a nonzero sized truth event collection, disabling truthMatching");
    doTruthMatching = false;
  }
  if(doTruthMatching) hardScatterEvent = mcEventCollection->at(0);

  if(p_SCT_RDO_cont.isValid()) {
    // loop over RDO container
    SCT_RDO_Container::const_iterator rdoCont_itr(p_SCT_RDO_cont->begin());
    const SCT_RDO_Container::const_iterator rdoCont_end(p_SCT_RDO_cont->end());

    for ( ; rdoCont_itr != rdoCont_end; ++rdoCont_itr ) {
      const SCT_RDO_Collection* p_SCT_RDO_coll(*rdoCont_itr);
      SCT_RDO_Collection::const_iterator rdo_itr(p_SCT_RDO_coll->begin());
      const SCT_RDO_Collection::const_iterator rdo_end(p_SCT_RDO_coll->end());

      for ( ; rdo_itr != rdo_end; ++rdo_itr ) {
        if(doTruthMatching){
          m_h_truthMatchedRDOs->Fill(1.5);
          bool findMatch = false;
          if(simDataMapSCT.isValid()){
            InDetSimDataCollection::const_iterator iter = (*simDataMapSCT).find((*rdo_itr)->identify());

            if ( iter != (*simDataMapSCT).end() ) {
              const InDetSimData& sdo = iter->second;
              const std::vector< InDetSimData::Deposit >& deposits = sdo.getdeposits();
              std::vector< InDetSimData::Deposit >::const_iterator nextdeposit = deposits.begin();
              std::vector< InDetSimData::Deposit >::const_iterator lastdeposit = deposits.end();
              for( ; nextdeposit!=lastdeposit; ++nextdeposit) {
	              const HepMcParticleLink& particleLink = nextdeposit->first;
                if(particleLink.isValid() && !findMatch){
                  HepMC::ConstGenParticlePtr genPart(particleLink.cptr());
                  if(genPart->parent_event() == hardScatterEvent) m_h_truthMatchedRDOs->Fill(3.5);
                  m_h_truthMatchedRDOs->Fill(2.5);
                  findMatch = true;
                }
              }
            }
          }
          if(!findMatch) m_h_truthMatchedRDOs->Fill(4.5);
        }
        const Identifier rdoID((*rdo_itr)->identify());
        const unsigned int rdoWord((*rdo_itr)->getWord());
        const int sctBrlEc(m_sctID->barrel_ec(rdoID));
        const int sctLayerDisk(m_sctID->layer_disk(rdoID));
        const int sctPhiMod(m_sctID->phi_module(rdoID));
        const int sctEtaMod(m_sctID->eta_module(rdoID));
        const int sctSide(m_sctID->side(rdoID));
        const int sctStrip(m_sctID->strip(rdoID));
	const int sctRow(m_sctID->row(rdoID));
        const int sctGroupSize((*rdo_itr)->getGroupSize());

        const unsigned long long rdoID_int = rdoID.get_compact();


	const InDetDD::SiDetectorElement *detEl = m_SCT_Manager->getDetectorElement(rdoID);

	if(!detEl) {
	  ATH_MSG_WARNING("No Element found for ID "<<m_sctID->show_to_string(rdoID)<<" - skipping!");
	  continue;
	}

	//Using lorentz-angle corrected version here will result in inconsistencies (ATLSWUPGR-103)
	//NB rationlization as here should probably also be done in master
	//https://gitlab.cern.ch/atlas/athena/-/merge_requests/33398

	Amg::Vector2D localPos = detEl->rawLocalPositionOfCell(rdoID);

	std::pair<Amg::Vector3D, Amg::Vector3D> endsOfStrip = detEl->endsOfStrip(localPos);

	if (m_doPosition) {

          m_globalX0->push_back(endsOfStrip.first.x());
          m_globalY0->push_back(endsOfStrip.first.y());
          m_globalZ0->push_back(endsOfStrip.first.z());

          m_globalX1->push_back(endsOfStrip.second.x());
          m_globalY1->push_back(endsOfStrip.second.y());
          m_globalZ1->push_back(endsOfStrip.second.z());

          m_localX->push_back(localPos.x());
          m_localY->push_back(localPos.y());
          m_localZ->push_back(0.0);

	}
        float stripradius0 = sqrt(endsOfStrip.first.x()*endsOfStrip.first.x()+endsOfStrip.first.y()*endsOfStrip.first.y());
        float stripradius1 = sqrt(endsOfStrip.second.x()*endsOfStrip.second.x()+endsOfStrip.second.y()*endsOfStrip.second.y());

        m_h_globalXY->Fill(endsOfStrip.first.x(),endsOfStrip.first.y());
        m_h_globalXY->Fill(endsOfStrip.second.x(),endsOfStrip.second.y());
        m_h_globalXY_shared->Fill(endsOfStrip.first.x(),endsOfStrip.first.y());
        m_h_globalXY_shared->Fill(endsOfStrip.second.x(),endsOfStrip.second.y());
        m_h_globalZR->Fill(endsOfStrip.first.z(),stripradius0);
        m_h_globalZR->Fill(endsOfStrip.second.z(),stripradius1);
        m_h_globalZR_shared->Fill(endsOfStrip.first.z(),stripradius0);
        m_h_globalZR_shared->Fill(endsOfStrip.second.z(),stripradius1);
        m_h_globalX->Fill(endsOfStrip.first.x());
        m_h_globalY->Fill(endsOfStrip.first.y());
        m_h_globalZ->Fill(endsOfStrip.first.z());
        m_h_globalX->Fill(endsOfStrip.second.x());
        m_h_globalY->Fill(endsOfStrip.second.y());
        m_h_globalZ->Fill(endsOfStrip.second.z());

        m_rdoID->push_back(rdoID_int);
        m_rdoWord->push_back(rdoWord);
        m_barrelEndcap->push_back(sctBrlEc);
        m_layerDisk->push_back(sctLayerDisk);
        m_phiModule->push_back(sctPhiMod);
        m_etaModule->push_back(sctEtaMod);
        m_side->push_back(sctSide);
	m_row->push_back(sctRow);
        m_strip->push_back(sctStrip);
        m_groupSize->push_back(sctGroupSize);

        m_h_rdoID->Fill(rdoID_int);
        m_h_rdoWord->Fill(rdoWord);
        m_h_barrelEndcap->Fill(sctBrlEc);
        m_h_layerDisk->Fill(sctLayerDisk);
        m_h_phiModule->Fill(sctPhiMod);
        m_h_etaModule->Fill(sctEtaMod);
        m_h_side->Fill(sctSide);
        m_h_strip->Fill(sctStrip);
	m_h_row->Fill(sctRow);
        m_h_groupSize->Fill(sctGroupSize);
        m_h_phi_v_eta->Fill(sctEtaMod, sctPhiMod);

        if (sctBrlEc == 0) {
          m_h_brlLayer->Fill(sctLayerDisk);
          m_h_brlPhiMod->Fill(sctPhiMod);
          m_h_brlEtaMod->Fill(sctEtaMod);
          m_h_brlSide->Fill(sctSide);
          m_h_brlStrip->Fill(sctStrip);
          m_h_brlGroupSize->Fill(sctGroupSize);
          m_h_brl_phi_v_eta->Fill(sctEtaMod, sctPhiMod);
          m_h_brl_strip_perLayer[sctLayerDisk]->Fill(sctStrip);
          m_h_brl_strip_perLayer[sctLayerDisk]->Fill(sctStrip);
        }
        else if (abs(sctBrlEc) == 2) {
          m_h_ecDisk->Fill(sctLayerDisk);
          m_h_ecPhiMod->Fill(sctPhiMod);
          m_h_ecEtaMod->Fill(sctEtaMod);
          m_h_ecSide->Fill(sctSide);
          m_h_ecStrip->Fill(sctStrip);
          m_h_ecGroupSize->Fill(sctGroupSize);
          m_h_ec_phi_v_eta->Fill(sctEtaMod, sctPhiMod);
          m_h_ec_strip_perLayer[sctLayerDisk]->Fill(sctStrip);
          m_h_ec_strip_perLayer[sctLayerDisk]->Fill(sctStrip);
        }
      }
    }
  }

  // SimData
  if(simDataMapSCT.isValid()) {
    // loop over SDO container
    InDetSimDataCollection::const_iterator sdo_itr(simDataMapSCT->begin());
    const InDetSimDataCollection::const_iterator sdo_end(simDataMapSCT->end());

    std::vector<int> barcode_vec;
    std::vector<int> eventIndex_vec;
    std::vector<float> charge_vec;
    for ( ; sdo_itr != sdo_end; ++sdo_itr ) {
      const Identifier sdoID((*sdo_itr).first);
      const InDetSimData& sdo((*sdo_itr).second);
      const unsigned long long sdoID_int = sdoID.get_compact();
      const int sdoWord(sdo.word());
      const int sctBrlEc_sdo(m_sctID->barrel_ec(sdoID));
      const int sctLayerDisk_sdo(m_sctID->layer_disk(sdoID));
      const int sctPhiMod_sdo(m_sctID->phi_module(sdoID));
      const int sctEtaMod_sdo(m_sctID->eta_module(sdoID));
      const int sctSide_sdo(m_sctID->side(sdoID));
      const int sctStrip_sdo(m_sctID->strip(sdoID));
      const int sctRow_sdo(m_sctID->row(sdoID));
      const bool noise(SCT_SimHelper::isNoise(sdo));
      const bool belowThresh(SCT_SimHelper::isBelowThreshold(sdo));
      const bool disabled(SCT_SimHelper::isDisabled(sdo));

      m_sdoID->push_back(sdoID_int);
      m_sdoWord->push_back(sdoWord);
      m_barrelEndcap_sdo->push_back(sctBrlEc_sdo);
      m_layerDisk_sdo->push_back(sctLayerDisk_sdo);
      m_phiModule_sdo->push_back(sctPhiMod_sdo);
      m_etaModule_sdo->push_back(sctEtaMod_sdo);
      m_side_sdo->push_back(sctSide_sdo);
      m_strip_sdo->push_back(sctStrip_sdo);
      m_row_sdo->push_back(sctRow_sdo);
      m_noise->push_back(noise);
      m_belowThresh->push_back(belowThresh);
      m_disabled->push_back(disabled);

      if (belowThresh) {
        if (sctBrlEc_sdo==0)
          m_h_belowThresh_brl->Fill(sctLayerDisk_sdo);
        else if (abs(sctBrlEc_sdo)==2)
          m_h_belowThresh_ec->Fill(sctLayerDisk_sdo);
      }

      if (disabled) {
        if (sctBrlEc_sdo==0)
          m_h_disabled_brl->Fill(sctLayerDisk_sdo);
        else if (abs(sctBrlEc_sdo)==2)
          m_h_disabled_ec->Fill(sctLayerDisk_sdo);
      }

      m_h_sdoID->Fill(sdoID_int);
      m_h_sdoWord->Fill(sdoWord);
      m_h_barrelEndcap_sdo->Fill(sctBrlEc_sdo);
      m_h_layerDisk_sdo->Fill(sctLayerDisk_sdo);
      m_h_phiModule_sdo->Fill(sctPhiMod_sdo);
      m_h_etaModule_sdo->Fill(sctEtaMod_sdo);
      m_h_side_sdo->Fill(sctSide_sdo);
      m_h_strip_sdo->Fill(sctStrip_sdo);
      m_h_row_sdo->Fill(sctRow_sdo);
      m_h_phi_v_eta_sdo->Fill(sctEtaMod_sdo, sctPhiMod_sdo);

      // loop over deposits
      // InDetSimData::Deposit typedef for std::pair<HepMCParticleLink,float>
      const std::vector<InDetSimData::Deposit>& deposits = sdo.getdeposits();
      std::vector<InDetSimData::Deposit>::const_iterator dep_itr(deposits.begin());
      const std::vector<InDetSimData::Deposit>::const_iterator dep_end(deposits.end());

      for ( ; dep_itr != dep_end; ++dep_itr ) {
        const HepMcParticleLink& particleLink = (*dep_itr).first;
        const int bar(particleLink.barcode());
        const int eventIx(particleLink.eventIndex());
        const float charge((*dep_itr).second);

        m_barcode->push_back(bar);
        m_eventIndex->push_back(eventIx);
        m_charge->push_back(charge);

        m_h_barcode->Fill(bar);
        m_h_eventIndex->Fill(eventIx);
        m_h_charge->Fill(charge);

        barcode_vec.push_back(bar);
        eventIndex_vec.push_back(eventIx);
        charge_vec.push_back(charge);
      }
      m_barcode_vec->push_back(barcode_vec);
      m_eventIndex_vec->push_back(eventIndex_vec);
      m_charge_vec->push_back(charge_vec);
      barcode_vec.clear();
      eventIndex_vec.clear();
      charge_vec.clear();
    }
  }

  if (m_tree) {
    m_tree->Fill();
  }

  return StatusCode::SUCCESS;
}

} // namespace ITk
