/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


#include "BCM_RDOAnalysis.h"

#include "TTree.h"
#include "TString.h"

#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>

BCM_RDOAnalysis::BCM_RDOAnalysis(const std::string& name, ISvcLocator *pSvcLocator)
  : AthAlgorithm(name, pSvcLocator)
  , m_word1(0)
  , m_word2(0)
  , m_chan(0)
  , m_pulse1Pos(0)
  , m_pulse1Width(0)
  , m_pulse2Pos(0)
  , m_pulse2Width(0)
  , m_LVL1A(0)
  , m_BCID(0)
  , m_LVL1ID(0)
  , m_err(0)
  , m_sdoID(0)
  , m_sdoWord(0)
  , m_barcode(0)
  , m_eventIndex(0)
  , m_charge(0)
  , m_barcode_vec(0)
  , m_eventIndex_vec(0)
  , m_charge_vec(0)

  , h_word1(0)
  , h_word2(0)
  , h_chan(0)
  , h_pulse1Pos(0)
  , h_pulse1Width(0)
  , h_pulse2Pos(0)
  , h_pulse2Width(0)
  , h_sdoID(0)
  , h_sdoWord(0)
  , h_barcode(0)
  , h_eventIndex(0)
  , h_charge(0)
    
  , m_tree(0)
  , m_ntupleFileName("/ntuples/file1")
  , m_ntupleDirName("/BCM_RDOAnalysis/")
  , m_ntupleTreeName("BCM_RDOAna")
  , m_path("/BCM_RDOAnalysis/")
  , m_thistSvc("THistSvc", name)
{
  declareProperty("NtupleFileName", m_ntupleFileName);
  declareProperty("NtupleDirectoryName", m_ntupleDirName);
  declareProperty("NtupleTreeName", m_ntupleTreeName);
  declareProperty("HistPath", m_path);
}

StatusCode BCM_RDOAnalysis::initialize() {
  ATH_MSG_DEBUG( "Initializing BCM_RDOAnalysis" );

  // Grab Ntuple and histgoramming service for tree
  CHECK(m_thistSvc.retrieve());

  m_tree = new TTree(TString(m_ntupleTreeName), "BCM_RDOAna");
  std::string fullNtupleName = "/" + m_ntupleFileName + "/" + m_ntupleDirName + "/" + m_ntupleTreeName;
  CHECK(m_thistSvc->regTree(fullNtupleName,m_tree));
  if (m_tree) {
    // BCM RDO
    m_tree->Branch("word1", &m_word1);
    m_tree->Branch("word2", &m_word2);
    m_tree->Branch("chan", &m_chan);
    m_tree->Branch("pulse1Pos", &m_pulse1Pos);
    m_tree->Branch("pulse1Width", &m_pulse1Width);
    m_tree->Branch("pulse2Pos", &m_pulse2Pos);
    m_tree->Branch("pulse2Width", &m_pulse2Width);
    m_tree->Branch("LVL1A", &m_LVL1A);
    m_tree->Branch("BCID", &m_BCID);
    m_tree->Branch("LVL1ID", &m_LVL1ID);
    m_tree->Branch("err", &m_err);
    // BCM SDO
    m_tree->Branch("sdoID", &m_sdoID);
    m_tree->Branch("sdoWord", &m_sdoWord);
    m_tree->Branch("barcode", &m_barcode);
    m_tree->Branch("eventIndex", &m_eventIndex);
    m_tree->Branch("charge", &m_charge);
    m_tree->Branch("barcode_vec", &m_barcode_vec);
    m_tree->Branch("eventIndex_vec", &m_eventIndex_vec);
    m_tree->Branch("charge_vec", &m_charge_vec);
  }
  else {
    ATH_MSG_ERROR( "No tree found" );
  }

  // HISTOGRAMS
  h_word1 = new TH1F("h_word1", "word 1", 100, 0, 4.5e7);
  h_word1->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_word1->GetName(), h_word1));

  h_word2 = new TH1F("h_word2", "word 2", 100, 0, 10);
  h_word2->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_word2->GetName(), h_word2));

  h_chan = new TH1F("h_chan", "channel ID", 100, 0, 20);
  h_chan->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_chan->GetName(), h_chan));

  h_pulse1Pos = new TH1F("h_pulse1Pos", "pulse 1 position", 100, 0, 70);
  h_pulse1Pos->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_pulse1Pos->GetName(), h_pulse1Pos));

  h_pulse1Width = new TH1F("h_pulse1Width", "pulse 1 width", 100, 0, 30);
  h_pulse1Width->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_pulse1Width->GetName(), h_pulse1Width));

  h_pulse2Pos = new TH1F("h_pulse2Pos", "pulse 2 position", 100, 0, 70);
  h_pulse2Pos->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_pulse2Pos->GetName(), h_pulse2Pos));

  h_pulse2Width = new TH1F("h_pulse2Width", "pulse 2 width", 100, 0, 30);
  h_pulse2Width->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_pulse2Width->GetName(), h_pulse2Width));

  h_sdoID = new TH1F("h_sdoID", "sdoID", 100, 0, 10);
  h_sdoID->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_sdoID->GetName(), h_sdoID));

  h_sdoWord = new TH1F("h_sdoWord", "sdoWord", 100, 0, 10);
  h_sdoWord->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_sdoWord->GetName(), h_sdoWord));

  h_barcode = new TH1F("h_barcode", "Barcode (SDO)", 100, 0, 2.2e9);
  h_barcode->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_barcode->GetName(), h_barcode));

  h_eventIndex = new TH1F("h_eventIndex", "Event index (SDO)", 100, 0, 1000);
  h_eventIndex->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_eventIndex->GetName(), h_eventIndex));

  h_charge = new TH1F("h_charge", "Charge (SDO)", 100, 0, 10);
  h_charge->StatOverflows();
  CHECK(m_thistSvc->regHist(m_path + h_charge->GetName(), h_charge));

  return StatusCode::SUCCESS;
  
}

StatusCode BCM_RDOAnalysis::execute() {
  ATH_MSG_DEBUG( "In BCM_RDOAnalysis::execute()" );

  m_word1->clear();
  m_word2->clear();
  m_chan->clear();
  m_pulse1Pos->clear();
  m_pulse1Width->clear();
  m_pulse2Pos->clear();
  m_pulse2Width->clear();
  m_LVL1A->clear();
  m_BCID->clear();
  m_LVL1ID->clear();
  m_err->clear();
  m_sdoID->clear();
  m_sdoWord->clear();
  m_barcode->clear();
  m_eventIndex->clear();
  m_charge->clear();
  m_barcode_vec->clear();
  m_eventIndex_vec->clear();
  m_charge_vec->clear();
  
  // RDO
  const BCM_RDO_Container* p_BCM_RDO_cont;
  if (evtStore()->retrieve(p_BCM_RDO_cont, "BCM_RDOs") == StatusCode::SUCCESS) {
    // loop over RDO container
    BCM_RDO_Container::const_iterator rdoCont_itr(p_BCM_RDO_cont->begin());
    const BCM_RDO_Container::const_iterator rdoCont_end(p_BCM_RDO_cont->end());

    for ( ; rdoCont_itr != rdoCont_end; ++rdoCont_itr ) {
      const DataVector<BCM_RawData>* p_BCM_RDO_coll(*rdoCont_itr);
      DataVector<BCM_RawData>::const_iterator rdo_itr(p_BCM_RDO_coll->begin());
      const DataVector<BCM_RawData>::const_iterator rdo_end(p_BCM_RDO_coll->end());

      for ( ; rdo_itr != rdo_end; ++rdo_itr ) {
	const int word1((*rdo_itr)->getWord1());
	const int word2((*rdo_itr)->getWord2());
	const int chan((*rdo_itr)->getChannel());
	const int pulse1Pos((*rdo_itr)->getPulse1Position());
	const int pulse1Width((*rdo_itr)->getPulse1Width());
	const int pulse2Pos((*rdo_itr)->getPulse2Position());
	const int pulse2Width((*rdo_itr)->getPulse2Width());
	const int LVL1A((*rdo_itr)->getLVL1A());
	const int BCID((*rdo_itr)->getBCID());
	const int LVL1ID((*rdo_itr)->getLVL1ID());
	const int err((*rdo_itr)->getError());
	
	m_word1->push_back(word1);
	m_word2->push_back(word2);
	m_chan->push_back(chan);
	m_pulse1Pos->push_back(pulse1Pos);
	m_pulse1Width->push_back(pulse1Width);
	m_pulse2Pos->push_back(pulse2Pos);
	m_pulse2Width->push_back(pulse2Width);
	m_LVL1A->push_back(LVL1A);
	m_BCID->push_back(BCID);
	m_LVL1ID->push_back(LVL1ID);
	m_err->push_back(err);

	h_word1->Fill(word1);
	h_word2->Fill(word2);
	h_chan->Fill(chan);
	h_pulse1Pos->Fill(pulse1Pos);
	h_pulse1Width->Fill(pulse1Width);
	h_pulse2Pos->Fill(pulse2Pos);
	h_pulse2Width->Fill(pulse2Width);
      }
    }
  }

  // SDO
  const InDetSimDataCollection* simDataMapBCM(NULL);
  if (evtStore()->retrieve(simDataMapBCM, "BCM_SDO_Map") == StatusCode::SUCCESS) {
    // loop over SDO container
    InDetSimDataCollection::const_iterator sdo_itr(simDataMapBCM->begin());
    const InDetSimDataCollection::const_iterator sdo_end(simDataMapBCM->end());

    std::vector<int> barcode_vec;
    std::vector<int> eventIndex_vec;
    std::vector<float> charge_vec;
    for ( ; sdo_itr != sdo_end; ++sdo_itr ) {
      const Identifier sdoID((*sdo_itr).first);
      const InDetSimData& sdo((*sdo_itr).second);
      
      const unsigned long long sdoID_int = sdoID.get_compact();
      const int sdoWord(sdo.word());
      
      m_sdoID->push_back(sdoID_int);
      m_sdoWord->push_back(sdoWord);

      h_sdoID->Fill(sdoID_int);
      h_sdoWord->Fill(sdoWord);

      // loop over deposits
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

	h_barcode->Fill(bar);
	h_eventIndex->Fill(eventIx);
	h_charge->Fill(charge);
	
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

StatusCode BCM_RDOAnalysis::finalize() {

  return StatusCode::SUCCESS;

}
