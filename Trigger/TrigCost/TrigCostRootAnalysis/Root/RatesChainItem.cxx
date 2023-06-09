// Dear emacs, this is -*- c++ -*-

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// -------------------------------------------------------------
//  author: Tim Martin <Tim.Martin@cern.ch>
// -------------------------------------------------------------

// STL include(s):
#include <bitset>

// Local include(s):
#include "../TrigCostRootAnalysis/RatesChainItem.h"
#include "../TrigCostRootAnalysis/CounterBaseRates.h"
#include "../TrigCostRootAnalysis/Config.h"
#include "../TrigCostRootAnalysis/TrigXMLService.h"

// ROOT includes
#include <TError.h>
#include <TMath.h>

namespace TrigCostRootAnalysis {

  UInt_t RatesChainItem::s_chainCount = 0; // Each ChainItem gets a sequential random seed

  /**
   * Construct new RatesChainItem with a given prescale.
   */
  RatesChainItem::RatesChainItem(std::string _name, Int_t _level, Double_t _PS, Double_t _PSExpress) :
    m_name(_name),
    m_level(_level),
    m_PS(_PS), // Integer prescale
    m_PSWeight(1./m_PS), // Reciprocal of the prescale - this is the basic weight quantity for this ChainItem
    m_PSReduced(1.),  
    m_PSReducedWeight(1.),
    m_PSExpress(_PSExpress),
    m_PSExpressWeight(1./_PSExpress),
    m_extraEfficiency(1.),
    m_R( ++s_chainCount ),
    m_ID( s_chainCount ),
    m_bunchGroupType(kBG_UNSET),
    m_lumiExtrapolationMap(nullptr),
    m_passRaw(kFALSE),
    m_passPS(kFALSE),
    m_inEvent(kFALSE),
    m_doDirectPS(kFALSE),
    m_matchRandomToOnline(kFALSE),
    m_advancedLumiScaling(kFALSE),
    m_iAmRandom(kFALSE),
    m_triggerLogic(nullptr),
    m_proxy(nullptr),
    m_bufferMu(0)
  {
    if (Config::config().debug()) {
      Info("RatesChainItem::RatesChainItem","New ChainItem:%s, Level:%i PS:%f", m_name.c_str(), m_level, m_PS);
    }

    m_doEBWeighting = Config::config().getInt(kDoEBWeighting); // Cache for speed
    m_doDirectPS = Config::config().getInt(kDirectlyApplyPrescales); // Cache for speed
    m_advancedLumiScaling = Config::config().getInt(kDoAdvancedLumiScaling);

    if (m_PS <= 0.) m_PSWeight = 0.;
    m_forcePassRaw = (Bool_t) Config::config().getInt(kRatesForcePass);
    m_matchRandomToOnline = (Bool_t) Config::config().getInt(kMatchL1RandomToOnline);
  }

  /**
   * Look at myself and classify from their name what type of BG I trigger on
   */
  void RatesChainItem::classifyBunchGroup() {

    RatesChainItem* _toCheck = this;
    if (m_level > 1 && m_lower.size() == 1) _toCheck = *(m_lower.begin());

    if (_toCheck->getName().find("_BPTX") != std::string::npos || _toCheck->getName().find("_BGRP") != std::string::npos) { // Ignore the beam pickup & specialist triggers
      m_bunchGroupType = kBG_NONE;
    } else if (_toCheck->getName().find("_FIRSTEMPTY") != std::string::npos) {
      m_bunchGroupType = kBG_FIRSTEMPTY;
    } else if (_toCheck->getName().find("_EMPTY") != std::string::npos) {
      m_bunchGroupType = kBG_EMPTY;
    } else if (_toCheck->getName().find("_UNPAIRED_ISO") != std::string::npos) {
      m_bunchGroupType = kBG_UNPAIRED_ISO;
    } else if (_toCheck->getName().find("_UNPAIRED_NONISO") != std::string::npos) {
      m_bunchGroupType = kBG_UNPAIRED_NONISO;
    } else if (_toCheck->getName().find("_ABORTGAPNOTCALIB") != std::string::npos) {
      m_bunchGroupType = kBG_ABORTGAPNOTCALIB;
    } else if (_toCheck->getName().find("_CALREQ") != std::string::npos) {
      m_bunchGroupType = kBG_CALREQ;
    } else {
      m_bunchGroupType = kBG_FILLED;
    }

    if (Config::config().debug()) {
      Info("RatesChainItem::classifyBunchGroup","Item %s classified as %s.",
        getName().c_str(), BunchGroupNameStr[m_bunchGroupType].c_str() );
    }

  }

  /**
   * Look at myself and classify if I am a random seeded L1 or HLT
   * If RANDOM and L1, then the rate is independent of lumi. This item gets a lumi extrap. factor of 1
   * If RANDOM and HLT, then the L1 rate is fixed, and the only lumi extrapolation comes from the increase in <mu>
   * If NOT RANDOM, then simple linear extrapolation holds.
   */
  void RatesChainItem::classifyLumiAndRandom() {

    classifyBunchGroup();

    RatesChainItem* _toCheck = this;
    if (m_level > 1 && m_lower.size() == 1) _toCheck = *(m_lower.begin());

    if ( _toCheck->getName().find("_RD0") != std::string::npos ||
         _toCheck->getName().find("_RD1") != std::string::npos ||
         _toCheck->getName().find("_RD2") != std::string::npos ||
         _toCheck->getName().find("_RD3") != std::string::npos ||
         _toCheck->getName().find("_L1RD0") != std::string::npos ||
         _toCheck->getName().find("_L1RD1") != std::string::npos ||
         _toCheck->getName().find("_L1RD2") != std::string::npos ||
         _toCheck->getName().find("_L1RD3") != std::string::npos ) {
      m_iAmRandom = kTRUE;
      if (Config::config().debug()) Info("RatesChainItem::classifyLumiAndRandom","Item %s classified as random", getName().c_str());
    } else {
      m_iAmRandom = kFALSE;
    }

    if (m_advancedLumiScaling == kFALSE) { // Just do linear extrapolation

      m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorLinear); // Linear (including deadtime)
      
    } else { // m_advancedLumiScaling == kTRUE. Chains can have different extrapolation  

      Bool_t _specialCase1 = kFALSE, _specialCase3 = kFALSE;

      if      (m_iAmRandom == kTRUE && m_level == 1) _specialCase1 = kTRUE;
      else if (m_iAmRandom == kTRUE && m_level >  1) _specialCase3 = kTRUE;

      if (checkPatternNoLumiWeight(getName()) || _specialCase1) { // SPECIAL CASE #1

        m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorDeadtimeOnly);
        Config::config().addVecEntry(kListOfNoLumiWeightChains, getName());

      } else if (checkPatternNoMuLumiWeight(getName())) { // SPECIAL CASE #2

        m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorBunchOnly);
        Config::config().addVecEntry(kListOfNoMuLumiWeightChains, getName());

      } else if (checkPatternNoBunchLumiWeight(getName()) || _specialCase3) { // SPECIAL CASE #3

        m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorMuOnly);
        Config::config().addVecEntry(kListOfNoBunchLumiWeightChains, getName());

      } else if (checkPatternExponentialWithMu(getName())) { // SPECIAL CASE #4

        if (m_level == 1) { // We allow for different slopes at L1 and HLT
          m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorExpoL1);
        } else {
          m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorExpoHLT);
        }
        Config::config().addVecEntry(kListOfExpoMuLumiWeightChains, getName());

      } else if (m_bunchGroupType == kBG_EMPTY || m_bunchGroupType == kBG_FIRSTEMPTY) {

        // For empty BG, we scale with the *inverse* number of bunches, i.e. more bunches is less empty hence less empty rate
        // Fist empty is a bit of a fudge - assumes that NBunch is proportional to NTrain
        // No deadtime correction here
        m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorEmpty);

      } else if (m_bunchGroupType == kBG_UNPAIRED_ISO || m_bunchGroupType == kBG_UNPAIRED_NONISO || 
                 m_bunchGroupType == kBG_ABORTGAPNOTCALIB || m_bunchGroupType == kBG_CALREQ) {

        // We have no idea how the UNPAIRED items will change in the prediction - for now, zero extrapolation
        // Abort gap and calreq are always the same so again, no extrap here
        m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorUnity);

      } else {

        // The vast majority of chains will get this. This already includes the deadtime weight
        m_lumiExtrapolationMap = &(TrigXMLService::trigXMLService().m_lumiScalingFactorLinear);

      }
    }
  }

  /**
   * @return What this item needs to be scaled by to extrapolate its lumi to the target
   * @see RatesChainItem::classifyRandom
   * @param _lb The Lumi Block to return the extrapolation weight for
   * @param _disableEventLumiExtrapolation if extrapolation is disabled e.g. for UpgradeRatesMonitor which does this via overlay
   */
  Double_t RatesChainItem::getLumiExtrapolationFactor(UInt_t _lb, Bool_t _disableEventLumiExtrapolation) {
    if (_disableEventLumiExtrapolation) return 1.;
    return m_lumiExtrapolationMap->at(_lb); //  Simple. One number per run | Advanced. Different strategies per chain, varying per lumi block
  }

  /**
   * User can supply additional scaling factors which will alter the effective efficiency of this chain and hence the rate
   * @param _extraEfficiency Additional scaling factor which will scale the rate when the item fires
   */
  void RatesChainItem::setExtraEfficiency(Double_t _extraEfficiency) {
    m_extraEfficiency *= _extraEfficiency;
  }

  /**
   * Equiv to reciprocal of @see RatesChainItem::setExtraEfficiency
   * @param _reductionFactor Scale rate down by this factor
   */
  void RatesChainItem::setRateReductionFactor(Double_t _reductionFactor) {
    m_extraEfficiency *= 1. / _reductionFactor;
  }


  /**
   * For HLT items, each seeding L1 item should be linked here by passing its pointer.
   * Note we do not own the lower chainItem
   * @param _lower The pointer to another RatesChainItem which seeds this instance.
   */
  void RatesChainItem::addLower(RatesChainItem* _lower) {
    m_lower.insert(_lower);
  }

  /**
   * For L1 items a link to any HLT chanins seeded should be added here
   * @param _lower The pointer to another RatesChainItem which is seeded by this instance.
   */
  void RatesChainItem::addUpper(RatesChainItem* _upper) {
    m_upper.insert(_upper);
  }

  /**
   * Registers that a rates counter makes use of this ChainItem. We can use this info to speed up
   * execution by only processing the counters which we need to.
   * Note we do not own the CounterRates object
   * @param _clinet The pointer to a CounterRates object which makes use of this ChainItem to calculate the event weight.
   */
  void RatesChainItem::addCounter(CounterBaseRates* _client) {
    m_clients.insert( _client );
  }

  /**
   * @return A const iterator to the start of this counter's set of seeding counters
   */
  ChainItemSetIt_t RatesChainItem::getLowerStart() {
    return m_lower.begin();
  }

  /**
   * @return A const iterator to the end of this counter's set of seeding counters
   */
  ChainItemSetIt_t RatesChainItem::getLowerEnd() {
    return m_lower.end();
  }

  /**
   * @return A reference to the set of lower, seeding, items of this item
   */
  ChainItemSet_t& RatesChainItem::getLower() {
    return m_lower;
  }

  /**
   * @return A const iterator to the start of this counter's set of seeded counters
   */
  ChainItemSetIt_t RatesChainItem::getUpperStart() {
    return m_upper.begin();
  }

  /**
   * @return A const iterator to the end of this counter's set of seeded counters
   */
  ChainItemSetIt_t RatesChainItem::getUpperEnd() {
    return m_upper.end();
  }

  /**
   * @return A reference to the set of seeded items of this item
   */
  ChainItemSet_t& RatesChainItem::getUpper() {
    return m_upper;
  }

  /**
   * @param _find A chain item pointer to find in this chain item's set of seeding triggers.
   * @return kTRUE if this ChainItem has the supplied ChainItem listed as one of its lower, seeding items.
   */
  Bool_t RatesChainItem::getLowerContains(RatesChainItem* _find) {
    return static_cast<Bool_t>( m_lower.count( _find ) );
  }

  /**
   * @param _set Reference to a set of chain item pointers to test against.
   * @return kTRUE if *all* of the ChainItems supplied in _set are also listed as lower items of this ChainItem
   */
  Bool_t RatesChainItem::getLowerContainsAll( std::set<RatesChainItem*>& _set ) {
    for (ChainItemSetIt_t _it = _set.begin(); _it != _set.end(); ++_it) { // Check we contain all these
      if (getLowerContains( (*_it) ) == kFALSE) return kFALSE;
    }
    return kTRUE;
  }

  /**
   * @param _find A chain item pointer to find in this chain item's set of seeded triggers.
   * @return kTRUE if this ChainItem has the supplied ChainItem listed as one of its upper, seeded items.
   */
  Bool_t RatesChainItem::getUpperContains(RatesChainItem* _find) {
    return static_cast<Bool_t>( m_upper.count( _find ) );
  }

  /**
   * @param _set Reference to a set of chain item pointers to test against.
   * @return kTRUE if *all* of the ChainItems supplied in _set are also listed as upper items of this ChainItem
   */
  Bool_t RatesChainItem::getUpperContainsAll( std::set<RatesChainItem*>& _set ) {
    for (ChainItemSetIt_t _it = _set.begin(); _it != _set.end(); ++_it) { // Check we contain all these
      if (getUpperContains( (*_it) ) == kFALSE) return kFALSE;
    }
    return kTRUE;
  }

  /**
   * @return The configured prescale value
   */
  Double_t RatesChainItem::getPS() {
    return m_PS;
  }

  /**
   * Sets a new prescale value
   */
  void RatesChainItem::setPS(Double_t _PS) {
    m_PS = _PS;
    m_PSWeight = 1./m_PS;
    if (m_PS <= 0.) m_PSWeight = 0.;
  }

  /**
   * Sets a reduced PS value. This is the component of the prescale which is not coherent with other chains in the CPS group
   */
  void RatesChainItem::setPSReduced(Double_t _PSReduced) {
    m_PSReduced = _PSReduced;
    m_PSReducedWeight = 1./m_PSReduced;
    if (m_PSReduced <= 0.) m_PSReducedWeight = 0.;
  }

  /**
   * @return The chain item's name
   */
  const std::string& RatesChainItem::getName() {
    return m_name;
  }

  /**
   * @return The chain sequential internal ID
   */
  UInt_t RatesChainItem::getID() {
    return m_ID;
  }

  /**
   * @param _passRaw If this ChainItem passed raw in this event.
   * @param _counterSet Set of counters we will process, add to it counters that I influence. This is pass-by-reference and is modified.
   */
  void RatesChainItem::beginEvent(Bool_t _passRaw,  CounterBaseRatesSet_t& _counterSet) {
    m_passRaw = _passRaw;
    m_inEvent = kTRUE;
    _counterSet.insert( m_clients.begin(), m_clients.end() );

    if (m_doDirectPS) newRandomPS(); //TODO - check this is only used for DirectPS application. Saves many calls to TRandom3

    // For random seeded triggers where the HLT was re-run, we need to check that we only run over unbiased events in the sample
    if (m_matchRandomToOnline == kTRUE && m_iAmRandom == kTRUE) {
      if ( Config::config().getInt(kCurrentEventWasRandomOnline) == kFALSE ) {
        m_passRaw = kFALSE;
        m_inEvent = kFALSE;
      }
    }
  }

  /**
   * @param _eventTOBs A TOBAccumulator of all TOBs in this event or pseudo-event (simulated high pileup TOB overlay).
   * Note - this function call requires a TriggerLogic pointer to be set, this logic will be used against the set of TOBs
   */
  void RatesChainItem::beginEvent(TOBAccumulator* _eventTOBs) {
    m_inEvent = kTRUE;
    static Bool_t _largeJetWindow = Config::config().getInt(kUpgradeJetLargeWindow);

    // For random seeded triggers where the HLT was re-run, we need to check that we only run over unbiased events in the sample
    if (m_matchRandomToOnline == kTRUE && m_iAmRandom == kTRUE) {
      if ( Config::config().getInt(kCurrentEventWasRandomOnline) == kFALSE ) {
        m_passRaw = kFALSE;
        m_inEvent = kFALSE;
        return;
      }
    }

    m_bufferMu = _eventTOBs->mu();

    // Loop over logic
    m_passRaw = kTRUE; // Assume we passed, see if we didn't
    for (const TriggerCondition& _condition : getTriggerLogic()->conditions()) {

      if (_condition.m_type == kMissingEnergyString) {
        if (_eventTOBs->METOverflow() == kFALSE && _eventTOBs->MET() <= _condition.m_thresh) {
          m_passRaw = kFALSE;
          break;
        } 
      } else if (_condition.m_type == kEnergyString) {
        if (_eventTOBs->TEOverflow() == kFALSE && _eventTOBs->TE() <= _condition.m_thresh) {
          m_passRaw = kFALSE;
          break; 
        }
      } else if (_condition.m_type == kHTString) {
        if (_eventTOBs->HT() <= _condition.m_thresh) {
          m_passRaw = kFALSE;
          break;
        } 
      } else if (_condition.m_type == kMHTString) {
        if (_eventTOBs->MHT() <= _condition.m_thresh) {
          m_passRaw = kFALSE;
          break; 
        }
      } else {  // For EM/JET/TAU/MU

        UInt_t _tobsPassingCondition = 0;
        for (const auto& _tob : _eventTOBs->TOBs() ) {
          if (_tob.m_type != _condition.m_type) continue; // Incorrect type (EM/TAU/MU etc.). Don't discriminate on this one
          Float_t _et = _tob.m_et;
          if (_tob.m_type == kJetString && _largeJetWindow == kTRUE) _et = _tob.m_etLarge;
          // Energy too low ?
          if (_tob.m_type == kMuonString) {
            if (_et < _condition.m_thresh) continue; // Muons are at set thresholds so should be <
          } else {
            if (_et <= _condition.m_thresh) continue; // From testing on jets, really does seem to be <=
          }
          if (TMath::Abs(_tob.m_eta) * 10 < _condition.m_min) continue; // eta too low
          if (TMath::Abs(_tob.m_eta) * 10 > _condition.m_max) continue; // eta too high
          if (_condition.m_iso != 0) { // Check isolation bits (if conditions require isolation)
            std::bitset<5> _tobIso = _tob.m_iso;
            std::bitset<5> _conditionIso = _condition.m_iso;
            Bool_t _pass = kTRUE;
            for (UInt_t _b = 0; _b < 5; ++_b) {
              if (_conditionIso.test(_b) == kTRUE && _tobIso.test(_b) == kFALSE) _pass = kFALSE;
            }
            if (_pass == kFALSE) continue; // A required isolation bit was not found
          }
          ++_tobsPassingCondition; // All requirements met
          // Histogram
          if (_tob.m_type == kJetString) m_bufferJetRoIEta.push_back(_tob.m_eta);
          else if (_tob.m_type == kMuonString) m_bufferMuRoIEta.push_back(_tob.m_eta);
          else if (_tob.m_type == kEmString) m_bufferEmRoIEta.push_back(_tob.m_eta);
          else if (_tob.m_type == kTauString) m_bufferTauRoIEta.push_back(_tob.m_eta);
          //if (_tobsPassingCondition == _condition.m_multi) break; // Do we have enough TOBs passing this condition? Bail out if so, don't need more
        }
        if (_tobsPassingCondition < _condition.m_multi) {
          m_passRaw = kFALSE; // A condition was not satisfied :( all must be satisfied. We cannot accept this event.
          break;
        }

      }
    }

    //Info("RatesChainItem::beginEvent","%s applying logic to %i TOBs (passed - %i) MET is %f TE is %f", getName().c_str(), _eventTOBs->TOBs().size(), (Int_t)m_passRaw, _eventTOBs->MET(), _eventTOBs->TE());
  }

  /**
   * Used in Upgrade Rates mode - we plot the eta distribution of the thresholds we pass and the multiplicity
   */
  void RatesChainItem::fillHistograms(DataStore& _dataStore, Float_t _weight, Float_t _bunchWeight) {
    for (const Float_t _value : m_bufferJetRoIEta) _dataStore.store(kVarJetEta, _value, _weight);
    for (const Float_t _value : m_bufferMuRoIEta)  _dataStore.store(kVarMuEta, _value, _weight);
    for (const Float_t _value : m_bufferEmRoIEta)  _dataStore.store(kVarEmEta, _value, _weight);
    for (const Float_t _value : m_bufferTauRoIEta) _dataStore.store(kVarTauEta, _value, _weight);

    if (m_bufferJetRoIEta.size() > 0) _dataStore.store(kVarJetNThresh, m_bufferJetRoIEta.size(), _weight);
    if (m_bufferMuRoIEta.size() > 0)  _dataStore.store(kVarMuNThresh, m_bufferMuRoIEta.size() , _weight);
    if (m_bufferEmRoIEta.size() > 0)  _dataStore.store(kVarEmNThresh, m_bufferEmRoIEta.size() , _weight);
    if (m_bufferTauRoIEta.size() > 0) _dataStore.store(kVarTauNThresh, m_bufferTauRoIEta.size(), _weight);

    _dataStore.store(kVarMu, m_bufferMu, _weight);
    _dataStore.store(kVarBunchWeight, _bunchWeight, _weight); // What part of the extrapolation was explicitly due to change in number of bunches
  }


  /**
   * Reset all flags to zero
   */
  void RatesChainItem::endEvent() {
    m_passRaw = kFALSE;
    m_passPS = kFALSE;
    m_inEvent = kFALSE;

    m_bufferJetRoIEta.clear();
    m_bufferMuRoIEta.clear();
    m_bufferEmRoIEta.clear();
    m_bufferTauRoIEta.clear();
  }

  /**
   * Update the random prescale to a new value
   */
  void RatesChainItem::newRandomPS() {
    if (m_PS <= 0.) {
      m_passPS = kFALSE;
    } else {
      m_passPS = (m_R.Rndm() < m_PSWeight);
    }
  }

  /**
   * @return If this chain item was executed this event, regardless of whether or not it passed PS or passed Raw.
   * For L1 items, this means we need to check the bunch group
   */
  Bool_t RatesChainItem::getInEvent() {
    if (m_level == 1) {
      // If this L1 item passed then it must have been in the event, we only store L1s which pass so return this.
      // Also, if doEBWeighting is false we do not know what the bunch group was online
      if (m_inEvent == kTRUE || m_doEBWeighting == kFALSE) return m_inEvent;
      // Otherwise we check the bunchgroup
      return (Config::config().getInt(kCurrentEventBunchGroupID) == m_bunchGroupType);
    } else { // HLT
      return m_inEvent;
    }
  }

  /**
   * @return If this chain item passed (raw) in this event.
   */
  Bool_t RatesChainItem::getPassRaw() {
    if (m_forcePassRaw == kTRUE) return kTRUE;
    return m_passRaw;
  }

  /**
   * @return If this chain item passed its prescale in this event.
   */
  Bool_t RatesChainItem::getPassPS() {
    return m_passPS;
  }

  /**
   * @return 1/Prescale weighting factor for this event. This is scaled by an optional user supplied extra efficiency factor which can modulate the rate
   */
  Double_t RatesChainItem::getPSWeight(Bool_t _includeExpress) {
    if (m_proxy != nullptr) return m_proxy->getLastWeight();
    if (_includeExpress == kTRUE) return m_PSWeight * m_PSExpressWeight * m_extraEfficiency;
    return m_PSWeight * m_extraEfficiency;
  }

  /**
   * The reduced prescacle is the chain prescale divided by the lowest prescale of any chain within
   * the same coherent prescale group.
   * @return 1/PrescaleReduced weighting factor for this event. This is scaled by an optional user supplied extra efficiency factor which can modulate the rate
   */
  Double_t RatesChainItem::getPSReducedWeight(Bool_t _includeExpress) {
    if (_includeExpress == kTRUE) return m_PSReducedWeight * m_PSExpressWeight * m_extraEfficiency;
    return m_PSReducedWeight * m_extraEfficiency;
  }

  /**
   * @return Zero if this chain did not pass raw, else returns 1/Prescale
   */
  Double_t RatesChainItem::getPassRawOverPS(Bool_t _includeExpress) {
    if (getPassRaw() == kFALSE) return 0.;
    return getPSWeight(_includeExpress);
  }

  /**
   * The reduced prescacle is the chain prescale divided by the lowest prescale of any chain within
   * the same coherent prescale group.
   * @return Zero if this chain did not pass raw, else returns 1/PrescaleReduced
   */
  Double_t RatesChainItem::getPassRawOverPSReduced(Bool_t _includeExpress) {
    if (getPassRaw() == kFALSE) return 0.;
    return getPSReducedWeight(_includeExpress);
  }

  /**
   * @return If this chain passed both Raw and PS
   */
  Bool_t RatesChainItem::getPassRawAndPS() {
    return (getPassRaw() && getPassPS());
  }

  /**
   * @param _tl Use a TriggerLogic to generate the pass/fail for this chain
   * Note this is required to use void beginEvent(TOBAccumulator* _eventTOBs, CounterBaseRatesSet_t& _counterSet)
   * RatesChainItem object does not own the trigger logic.
   */
  void RatesChainItem::setTriggerLogic(TriggerLogic* _tl) {
    m_triggerLogic = _tl;
  }

  /**
   * @return Pointer to any registered TriggerLogic item
   */
  TriggerLogic* RatesChainItem::getTriggerLogic() {
    return m_triggerLogic;
  }

} // namespace TrigCostRootAnalysis


