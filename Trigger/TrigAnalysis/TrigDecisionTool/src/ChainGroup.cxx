/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**********************************************************************************
 * @Project: TrigDecisionTool
 * @Package: TrigDecisionTool
 * @Class  : ChainGroup
 *
 * @brief simple container to hold trigger chains
 *
 * @author Michael Begel  <michael.begel@cern.ch> - Brookhaven National Laboratory
 * @author Alexander Mann  <mann@cern.ch> - University of Goettingen
 *
 ***********************************************************************************/
#include <limits>
#include "boost/regex.hpp"
#include "boost/foreach.hpp"

#include "TrigConfHLTData/HLTChain.h"
#include "TrigConfHLTData/HLTStreamTag.h"
#include "TrigConfHLTData/HLTSignature.h"
#include "TrigConfHLTData/HLTTriggerElement.h"
#include "TrigConfHLTData/HLTUtils.h"
#include "TrigConfL1Data/TriggerItem.h"
#include "TrigConfL1Data/TriggerItemNode.h"
#include "TrigSteeringEvent/Chain.h"
#include "TrigSteeringEvent/Lvl1Result.h"
#include "TrigNavigation/ComboIterator.h"

#include "TrigDecisionTool/ChainGroup.h"
#include "TrigDecisionTool/TDTUtilities.h"
#include "TrigDecisionTool/Logger.h"

#define MLOG(x)   if (log().level()<=MSG::x) log() << MSG::x

using namespace std;

Trig::ChainGroup::ChainGroup(const std::vector< std::string >& triggerNames,
                             const Trig::CacheGlobalMemory&  parent) 
  : Logger("ChainGroup"),
    m_patterns(triggerNames),
    m_cgm(parent),
    m_prescale(0.)
{}

const Trig::ChainGroup& Trig::ChainGroup::operator+(const Trig::ChainGroup& rhs) {
  std::vector< std::string > v;
  v.resize(patterns().size()+rhs.patterns().size());
  merge(patterns().begin(),patterns().end(),
	(rhs.patterns()).begin(),(rhs.patterns()).end(),
	v.begin());
  return *(const_cast<Trig::CacheGlobalMemory*>(cgm())->createChainGroup(v));
}


bool Trig::ChainGroup::operator==(const Trig::ChainGroup& rhs) {
  std::vector< std::string > v1=Trig::keyWrap(names());
  std::vector< std::string > v2=Trig::keyWrap(rhs.names());
  if (v1==v2) return true;
  return false;
}

bool Trig::ChainGroup::operator!=(const Trig::ChainGroup& rhs) {
  std::vector< std::string > v1=Trig::keyWrap(names());
  std::vector< std::string > v2=Trig::keyWrap(rhs.names());
  if (v1!=v2) return true;
  return false;
}

void Trig::ChainGroup::addAlias(const std::string& alias) {
  const_cast<Trig::CacheGlobalMemory*>(cgm(true))->createChainGroup(patterns(),alias);
}

const Trig::CacheGlobalMemory* Trig::ChainGroup::cgm(bool onlyConfig) const {
  if ( !onlyConfig )
    const_cast<Trig::CacheGlobalMemory*>(&m_cgm)->assert_decision(); 
  
  return &m_cgm;
}

bool Trig::ChainGroup::HLTResult(const std::string& chain, unsigned int condition) const {
  bool chainRESULT = false;
  if (chain=="") return chainRESULT;
  const HLT::Chain* fchain=cgm()->chain(chain);
  if (fchain==0) return chainRESULT;
  
  
  bool RAW =         fchain->chainPassedRaw();
  bool PASSTHROUGH = fchain->isPassedThrough();
  bool PRESCALED   = fchain->isPrescaled();
  bool RESURRECTED = fchain->isResurrected();
  
  // Resurrection overwrites the value in RAW but sets the RESURRECTED flag
  // we should therefore fix RAW appropriately
  if (~condition & TrigDefs::allowResurrectedDecision) {
    if (RESURRECTED) {
      RAW=false;
    }
  }
  //
  // Do we accept the result?
  //
  if (condition & TrigDefs::passedThrough) {
    if (PASSTHROUGH) {chainRESULT=true;}
  }
  if (condition & TrigDefs::requireDecision) {
    if (RAW && !PRESCALED) {chainRESULT=true;}
    if ( condition & TrigDefs::allowResurrectedDecision ) { // prescaling does not matter for RR (it runs in fact because of that
      if (RAW) {chainRESULT=true;}      
    }

  }
  // respects resurrection -- is this the appropriate behavior???
  if (condition & TrigDefs::eventAccepted) {
    if ( (RAW  && !PRESCALED) ||  PASSTHROUGH) {chainRESULT=true;}
  }
  MLOG(DEBUG) << "ChainGroup::HLTResult Counter = " << std::setw(4) << fchain->getChainCounter()
	      << " name = "  << fchain->getChainName()
	      << " level = " << fchain->getConfigChain()->level()
	      << " success (raw) = " << fchain->chainPassedRaw()
	      << " pass-through = " << fchain->isPassedThrough()
	      << " prescaled = " << fchain->isPrescaled()
	      << " rerun = " << fchain->isResurrected()
	      << " lastActiveStep = " << fchain->getChainStep()
	      << " name = " << std::setw(35) << fchain->getChainName() 
	      << " result = " << chainRESULT
	      << endreq; 

  return chainRESULT;
}

// this logic fails for passthrough especially with enforceLogicalFlow!!!!
bool Trig::ChainGroup::L1Result(const std::string& item, unsigned int condition) const {
  bool r = false;
  if (item=="") return r;
  if (item.find(',')!=std::string::npos) {
    std::vector< std::string > items = convertStringToVector(item);
    std::vector< std::string >::iterator itit = items.begin();
    for(;itit != items.end(); ++itit)
      if(L1Result(*itit,condition)) return true;
    return false;
  }
  const LVL1CTP::Lvl1Item* fitem=cgm()->item(item);
  if (fitem==0) return r;
  MLOG(DEBUG) << " success (raw) = " << fitem->isPassedBeforePrescale()
	<< " prescaled = " << fitem->isPrescaled()
	<< " vetoed = " << fitem->isVeto()
	<< " name = " << std::setw(35) << fitem->name() 
	<< endreq; 

  r = fitem->isPassedAfterVeto();

  if (condition & TrigDefs::allowResurrectedDecision)
    r = fitem->isPassedBeforePrescale();

  //if (condition & TrigDefs::eventAccepted) 
  //  r = fitem->isPassedAfterVeto();
  //else
  //  r = fitem->isPassedBeforePrescale();
  return r;
}


std::string Trig::ChainGroup::getLowerName(const std::string& name) const {
  if ( name == "" )
    return "";
  const TrigConf::HLTChain* cchain = cgm(true)->config_chain(name);
  if (cchain==0){
    log() << MSG::WARNING
	  << " Lower chain name used by:  " << name << " is not in the configuration " << endreq;
    return "BAD NAME";
  }
  return cchain->lower_chain_name();
}

  

bool Trig::ChainGroup::isPassed(unsigned int condition) const 
{

  if (msgLvl(MSG::DEBUG)) {
    log() << MSG::DEBUG << " Got CG to work with " << patterns() << endreq;
    log() << MSG::DEBUG << " Got CG to work with " << names() << endreq;
  }

  ChainGroup::const_conf_chain_iterator chIt;

  bool RESULT = false;
  bool chainRESULT;
  std::string nexttwo;

  for ( chIt = conf_chain_begin(); chIt != conf_chain_end(); ++chIt) {
    chainRESULT=false;
    chainRESULT=HLTResult((*chIt)->chain_name(),condition);
    if (chainRESULT && (condition & TrigDefs::enforceLogicalFlow)) {
      // enforceLogicalFlow
      if ((*chIt)->level()=="EF") {
        nexttwo=getLowerName((*chIt)->chain_name());
        chainRESULT = chainRESULT && HLTResult(nexttwo,condition);
        chainRESULT = chainRESULT && L1Result(getLowerName(nexttwo),condition);

      } else if ((*chIt)->level()=="L2") {
        chainRESULT = chainRESULT && L1Result(getLowerName((*chIt)->chain_name()),condition);      

      } else if ((*chIt)->level()=="HLT"){
	    chainRESULT = chainRESULT && L1Result(getLowerName((*chIt)->chain_name()),condition);
      }
      
    }
    RESULT = RESULT || chainRESULT;
  }

  ChainGroup::const_conf_item_iterator iIt;
  for ( iIt = conf_item_begin(); iIt != conf_item_end(); ++iIt) {
    RESULT = RESULT || L1Result((*iIt)->name(),condition);
  }
  return RESULT;
}

unsigned int Trig::ChainGroup::HLTBits(const std::string& chain, const std::string& level) const {
  unsigned int chainRESULT = 0;
  if (chain=="") return chainRESULT;
  const HLT::Chain* fchain = cgm()->chain(chain);
  if (fchain==0) return chainRESULT;
  if (level=="L2") {
    if (fchain->chainPassedRaw())  chainRESULT = chainRESULT | TrigDefs::L2_passedRaw;
    if (fchain->isPassedThrough()) chainRESULT = chainRESULT | TrigDefs::L2_passThrough;
    if (fchain->isPrescaled())     chainRESULT = chainRESULT | TrigDefs::L2_prescaled;
    if (fchain->isResurrected())   chainRESULT = chainRESULT | TrigDefs::L2_resurrected;
  } else {//L2EF merged use same EF bits
    if (fchain->chainPassedRaw())  chainRESULT = chainRESULT | TrigDefs::EF_passedRaw;
    if (fchain->isPassedThrough()) chainRESULT = chainRESULT | TrigDefs::EF_passThrough;
    if (fchain->isPrescaled())     chainRESULT = chainRESULT | TrigDefs::EF_prescaled;
    if (fchain->isResurrected())   chainRESULT = chainRESULT | TrigDefs::EF_resurrected;
  }
  return chainRESULT;
}

unsigned int Trig::ChainGroup::L1Bits(const std::string& item) const {
  unsigned int r = 0;
  if (item=="") return r;
  if (item.find(',')!=std::string::npos) {
    std::vector< std::string > items = convertStringToVector(item);
    std::vector< std::string >::iterator itit = items.begin();
    for(;itit != items.end(); ++itit)
      r |= L1Bits(*itit);
    return r;
  }
  const LVL1CTP::Lvl1Item* fitem = cgm()->item(item);
  if (fitem==0) return r;
  if (fitem->isPassedBeforePrescale()) r = r | TrigDefs::L1_isPassedBeforePrescale;
  if (fitem->isPassedAfterPrescale())  r = r | TrigDefs::L1_isPassedAfterPrescale;
  if (fitem->isPassedAfterVeto())      r = r | TrigDefs::L1_isPassedAfterVeto;
  return r;
}

unsigned int Trig::ChainGroup::isPassedBits() const 
{

  ChainGroup::const_conf_chain_iterator chIt;

  unsigned int RESULT = 0;
  //  unsigned int chainRESULT=0;
  std::string nexttwo;

  for ( chIt = conf_chain_begin(); chIt != conf_chain_end(); ++chIt) {

    RESULT = RESULT | HLTBits((*chIt)->chain_name(),(*chIt)->level());

    if ((*chIt)->level()=="EF") {
      nexttwo=getLowerName((*chIt)->chain_name());
      RESULT = RESULT | HLTBits(nexttwo,"L2");
      RESULT = RESULT | L1Bits(getLowerName(nexttwo));
      
    } else if ((*chIt)->level()=="L2") {
      RESULT = RESULT | L1Bits(getLowerName((*chIt)->chain_name()));
      
    } else if ((*chIt)->level()=="HLT") {
      RESULT = RESULT | L1Bits(getLowerName((*chIt)->chain_name()));
    }
    RESULT = RESULT | RESULT;
    //    cout << "After looking at :" << (*chIt)->chain_name() 
    //	 << " " << std::hex << RESULT << endl;
  }


  ChainGroup::const_conf_item_iterator iIt;
  for ( iIt = conf_item_begin(); iIt != conf_item_end(); ++iIt) {
    RESULT = RESULT | L1Bits((*iIt)->name());
    //    cout << "After looking at :" << (*iIt)->name() << std::hex << RESULT << endl;
  }
  return RESULT;
}


HLT::ErrorCode Trig::ChainGroup::error() const {
  HLT::ErrorCode errorCode = HLT::OK;  
  ChainGroup::const_conf_chain_iterator chIt;
  for ( chIt = conf_chain_begin(); chIt != conf_chain_end(); ++chIt) {
    const HLT::Chain* fchain = cgm()->chain((*chIt)->chain_name());
    if (fchain==0) continue; 
    HLT::ErrorCode ec = fchain->getErrorCode();
    errorCode = errorCode > ec ? errorCode : ec;
  }
  return errorCode;
}

float Trig::ChainGroup::HLTPrescale(const std::string& chain, unsigned int /*condition*/) const {
  if (chain=="") return 0.;

  const TrigConf::HLTChain* fchain=cgm(true)->config_chain(chain);
  if (fchain==0) { // this is error condition, we always need configuration of the chains in the chaon group!
    log() << MSG::WARNING
          << "Configuration for the chain: " << chain << " not known" << endreq;
    return std::numeric_limits<float>::quiet_NaN();
  }
  float chainRESULT = fchain->prescale();
  
  //  log() << MSG::VERBOSE
  //  	<< "Configuration for the chain: " << chain << " for prescale: " << chainRESULT << endreq;

  /*
  // Resurrection makes this unprescaled!
  if (fchain->isResurrected() && condition & TrigDefs::allowResurrectedDecision) {
  const TrigConf::HLTChain* cchain=fchain->getConfigChain();
  if (cchain) {
  chainRESULT=cchain->rerun_prescale();
  }
  }
  */
  if (chainRESULT < 1)
    chainRESULT = 0.;
    
  return chainRESULT;
}


int Trig::ChainGroup::L1Prescale(const std::string& item, unsigned int /*condition*/) const {
  if (item=="") return 0;

  if(item.find(',')==std::string::npos) {
    const  TrigConf::TriggerItem* fitem=cgm(true)->config_item(item);
    if (fitem==0) {
      log() << MSG::WARNING
            << "Configuration for the item: " << item << " not known" << endreq;
      return std::numeric_limits<int>::quiet_NaN();
    }
    // now we can;t access the prescale value because this information doe not come togehther as in HLT
    // we need to go to the cache of L1 items and get it from there  
    int ctpid = fitem->ctpId();
    int itemprescale = cgm(true)->item_prescale(ctpid);
    if ( itemprescale < 1)
      itemprescale = 0;
    return itemprescale;
  } else {
    int minprescale=0;
    std::vector< std::string > items = convertStringToVector(item);
    std::vector< std::string >::iterator itit = items.begin();
    for(;itit != items.end(); ++itit) {
      const  TrigConf::TriggerItem* fitem=cgm(true)->config_item(*itit);
      if (fitem==0) {
        log() << MSG::WARNING << "Configuration for the item: " << *itit << " not known" << endreq;
        return std::numeric_limits<int>::quiet_NaN();
      }
      int ctpid = fitem->ctpId();
      int itemprescale = cgm(true)->item_prescale(ctpid);
      if ( itemprescale < 1)
	itemprescale = 0;
      minprescale = (minprescale&&(minprescale<itemprescale)?minprescale:itemprescale); // takes min, except the first time
    }
    return minprescale;
  }
}

float Trig::ChainGroup::getPrescale(unsigned int condition) const {
  if ( condition != TrigDefs::Physics )
    return 0.0;
  return m_prescale;
}

float Trig::ChainGroup::calculatePrescale(unsigned int condition)
{
  bool singleTrigger = (m_confChains.size()+m_confItems.size()==1);

  ChainGroup::const_conf_chain_iterator chIt;
  for ( chIt = conf_chain_begin(); chIt != conf_chain_end(); ++chIt) {

    const std::string & hltChainName = (*chIt)->chain_name();
    float chainRESULT = HLTPrescale(hltChainName,condition);

    if (condition & TrigDefs::enforceLogicalFlow) {
      // enforceLogicalFlow
      if ((*chIt)->level()=="EF") {
        std::string hltChainNameL2 = getLowerName(hltChainName);
        std::string l1ItemName     = getLowerName(hltChainNameL2);
        chainRESULT *= HLTPrescale(hltChainNameL2,condition);
        chainRESULT *= L1Prescale(l1ItemName,condition);
        if(l1ItemName.find(',')!=std::string::npos) singleTrigger=false;

      } else if ((*chIt)->level()=="L2") {
        std::string l1ItemName       = getLowerName(hltChainName);
        chainRESULT *= L1Prescale(l1ItemName,condition);
        if(l1ItemName.find(',')!=std::string::npos) singleTrigger=false;

	//todo: not clear how to handle prescales in the merged system??
      } else if ((*chIt)->level()=="HLT") {
        std::string l1ItemName       = getLowerName(hltChainName);
        chainRESULT *= L1Prescale(l1ItemName,condition);
        if(l1ItemName.find(',')!=std::string::npos) singleTrigger=false;

      }

    }

    if (singleTrigger) return chainRESULT;  // for a single trigger we are done

    bool UNPRESCALED = (fabs(chainRESULT-1.0)<1e-5);

    if (UNPRESCALED) return 1.0; // any unprescaled trigger and we are done too
  }


  ChainGroup::const_conf_item_iterator iIt;
  for ( iIt = conf_item_begin(); iIt != conf_item_end(); ++iIt) {

    const std::string & l1ItemName = (*iIt)->name();
    int itemRESULT = L1Prescale(l1ItemName,condition);
    if(l1ItemName.find(',')!=std::string::npos) singleTrigger=false;

    if (singleTrigger) return itemRESULT; // for a single trigger we are done

    bool UNPRESCALED = (itemRESULT==1);

    if (UNPRESCALED) return 1.0; // any unprescaled trigger and we are done too
  }

  return 0.0; // multiple triggers and all are prescaled
}


std::vector< std::string > Trig::ChainGroup::getListOfTriggers() const {
  std::vector< std::string > v;
  v.assign(m_names.begin(),m_names.end());
  return v;
}


std::vector< std::string > Trig::ChainGroup::getListOfStreams() const {
   std::set< std::string > s;
   std::set<const TrigConf::HLTChain*>::const_iterator foo;
   for (foo=conf_chain_begin(); foo != conf_chain_end(); ++foo) {
      const std::vector<TrigConf::HLTStreamTag*>& streamTagList = (*foo)->streams();
      std::vector<TrigConf::HLTStreamTag*>::const_iterator bar;
      for (bar=streamTagList.begin(); bar != streamTagList.end(); ++bar) {
         s.insert((*bar)->stream());
      }
   }
   std::vector< std::string > v;
   v.assign(s.begin(),s.end());
   return v;
}

//
// Groups
//

vector<string>
Trig::ChainGroup::getListOfGroups() const {

   vector< string > v;

   BOOST_FOREACH( const TrigConf::HLTChain* ch, m_confChains )
      v.assign( ch->groups().begin(), ch->groups().end() );

   return v;
}

//
// Signatures
//

std::vector< std::string > Trig::ChainGroup::getListOfSignatures() const {
  std::set< std::string > s;
  std::set<const TrigConf::HLTChain*>::const_iterator foo;
  // HLT Chain
  for (foo=conf_chain_begin(); foo != conf_chain_end(); ++foo) {
    const std::vector<TrigConf::HLTSignature*>& signatureList = (*foo)->signatureList();
    std::vector<TrigConf::HLTSignature*>::const_iterator bar;
    for (bar=signatureList.begin(); bar != signatureList.end(); ++bar) {
      s.insert((*bar)->label());
    }
  }
  std::vector< std::string > v;
  v.assign(s.begin(),s.end());
  return v;
}


//
// Return level 1 thresholds
//

vector<string>
Trig::ChainGroup::getListOfThresholds() const {

   set<string> s;   // using a set makes the items in the result vector unique
   std::stack<TrigConf::TriggerItemNode*> nodes;
   TrigConf::TriggerItemNode* node;

   BOOST_FOREACH( const TrigConf::TriggerItem* item, m_confItems ) {
      nodes.push( item->topNode() );
      while (!nodes.empty()) {
         node = nodes.top(); nodes.pop();
         if (node == NULL)
            continue;
         if (node->isThreshold()) {
            if (node->triggerThreshold()) {
               // available if thresholds have been read in
               if (!node->triggerThreshold()->name().empty())
                  s.insert(node->triggerThreshold()->name());
            } else if (!node->thresholdName().empty()) {
               // fall back solution
               s.insert(node->thresholdName());
            }
         } else {
            BOOST_FOREACH(TrigConf::TriggerItemNode* childnode, node->children()) {
               nodes.push(childnode);
            }
         }
      }
      // I am not using (*it)->topNode()->getAllThresholds() here, because it returns nothing when only ItemDef (and not the thresholds themselves) are defined
   }

   vector<string> v;
   v.assign(s.begin(),s.end());
   return v;

}


//
// Trigger Elements
//

std::vector< std::vector< std::string > > Trig::ChainGroup::getListOfTriggerElements() const {
  std::set< std::vector< std::string > > s;

  std::vector< std::string > t;

  std::set<const TrigConf::HLTChain*>::const_iterator foo;
  for (foo=conf_chain_begin(); foo != conf_chain_end(); ++foo) {
    const std::vector<TrigConf::HLTSignature*>& signatureList = (*foo)->signatureList();
    std::vector<TrigConf::HLTSignature*>::const_iterator bar;
    for (bar=signatureList.begin(); bar != signatureList.end(); ++bar) {
      t.clear();
      const std::vector<TrigConf::HLTTriggerElement*>& outputTEs = (*bar)->outputTEs();
      std::vector<TrigConf::HLTTriggerElement*>::const_iterator foobar;
      for (foobar=outputTEs.begin(); foobar != outputTEs.end(); ++foobar) {
	t.push_back( (*foobar)->name());
      }
      s.insert(t);
    }
  }
  std::vector< std::vector< std::string > > v;
  v.assign(s.begin(),s.end());
  return v;
}


//
// get vector of vector with all Trigger Elements
//

std::vector< std::vector< TrigConf::HLTTriggerElement* > > Trig::ChainGroup::getHLTTriggerElements() const {
  std::set< std::vector< TrigConf::HLTTriggerElement* > > s;

  std::vector< TrigConf::HLTTriggerElement* > t;

  std::set<const TrigConf::HLTChain*>::const_iterator foo;
  for (foo=conf_chain_begin(); foo != conf_chain_end(); ++foo) {
    const std::vector<TrigConf::HLTSignature*>& signatureList = (*foo)->signatureList();
    std::vector<TrigConf::HLTSignature*>::const_iterator bar;
    for (bar=signatureList.begin(); bar != signatureList.end(); ++bar) {
      t.clear();
      const std::vector<TrigConf::HLTTriggerElement*>& outputTEs = (*bar)->outputTEs();
      std::vector<TrigConf::HLTTriggerElement*>::const_iterator foobar;
      for (foobar=outputTEs.begin(); foobar != outputTEs.end(); ++foobar) {
	t.push_back( (*foobar) );
      }
      s.insert(t);
    }
  }
  std::vector< std::vector< TrigConf::HLTTriggerElement* > > v;
  v.assign(s.begin(),s.end());
  return v;
}


Trig::ChainGroup::~ChainGroup() {}

void
Trig::ChainGroup::update(const TrigConf::HLTChainList* confChains,
                         const TrigConf::ItemContainer* confItems) {

   m_confChains.clear();
   m_confItems.clear();
   m_names.clear();


   // protect against genConf failure
   if (!(confChains && confItems) ) return;

   for(std::vector< std::string >::const_iterator it = m_patterns.begin();
       it != m_patterns.end(); it++) {
      // find chains matching pattern     
      boost::regex compiled(*it);
      boost::cmatch what;

      BOOST_FOREACH(TrigConf::HLTChain* ch, *confChains) {
         if ( boost::regex_match(ch->chain_name().c_str(), what, compiled) ) {
            m_confChains.insert(ch);
            m_names.push_back(ch->chain_name());
         }
      }

      BOOST_FOREACH(TrigConf::TriggerItem* item, *confItems) {
         if ( boost::regex_match( item->name().c_str(), what, compiled) ) {
            m_confItems.insert(item);
            m_names.push_back(item->name());
         }
      }
   }
   m_prescale = calculatePrescale(TrigDefs::Physics);
}

/////////////////////////////////////////////////////////////////////////////
// features

namespace ChainGroup_impl {
   using namespace HLT;
   using namespace Trig;

   bool allActive(const std::vector<TriggerElement*>& tes) {
      BOOST_FOREACH(TriggerElement* te, tes){
         if (te->getActiveState() == false)
            return false;
      }
      return true;
   }


  
   void collectCombinations( const TrigConf::HLTChain* conf, const CacheGlobalMemory* cgm, FeatureContainer& fc, unsigned int condition ){
      // go over the steps of the chain and collecte TEs combinations for each of the chain step (signature)
      bool last_step=true;
      const TrigConf::HLTSignature* previous_sig(0);
      BOOST_REVERSE_FOREACH(const TrigConf::HLTSignature* sig, conf->signatureList()) {
         // chain without signatures
         if (!sig) break;

         // the signatures size changes from step to step, this needs custom treatement and the iteration eeds to be stoped here
         if ( previous_sig && previous_sig->outputTEs().size() != sig->outputTEs().size() )
            break;
         previous_sig = sig;

         std::vector<std::vector<HLT::TriggerElement*> > tes(sig->outputTEs().size()); // preallocate     
         size_t idx = 0;
         BOOST_FOREACH(const TrigConf::HLTTriggerElement* confte, sig->outputTEs() ) {

            // here the condition enters; if we take only accepted objects we need to pick only the active TEs
            cgm->navigation()->getAllOfType(confte->id(), tes[idx], condition & TrigDefs::Physics ); 	
            idx++;
         }
         HLT::ComboIterator combination(tes, cgm->navigation());

         // build the combinations, sometimes a huge list
         while (combination.isValid()) {
            // very very tricky, if all TEs in the combination are active then it means they were already picked up by previous combinations
            // but we can not do this for the last chain step, (we woudl be unable to pick objects whcih made trigger passing)
            //std::cerr << "emitting new combination last_step: "<< last_step << std::endl;	
            if (!allActive(*combination) || last_step) {
               fc.addWithChecking(Combination(*combination, cgm));
            }
            ++combination;
         }

         if ( condition & TrigDefs::Physics ) // here again athe condition calls for the loop breaking
            break;
         else
            last_step = false; // and go deeper
      }
   }
}// eof namespace


const Trig::FeatureContainer 
Trig::ChainGroup::features(unsigned int condition) const {
   using namespace ChainGroup_impl;
   FeatureContainer f(cgm());

   // this loop only applies to L2 and EF chain groups
   std::set<const TrigConf::HLTChain*>::const_iterator chIt;
   for (chIt=conf_chain_begin(); chIt != conf_chain_end(); ++chIt) { 
      const HLT::Chain* fchain = cgm()->chain(**chIt);
      if (fchain) {
         collectCombinations(fchain->getConfigChain(), cgm(), f, condition);
      }
   }


   // this part only applies to L1 chain groups
   std::vector< std::vector< HLT::TriggerElement*> > tes;
   std::vector< std::vector< HLT::TriggerElement*> >::iterator tesit;
   const_conf_item_iterator it; // iterator over configuration items, *it is TrigConf::TriggerItem*

   BOOST_FOREACH(const TrigConf::TriggerItem* item, m_confItems) {

      std::set< std::string > threshold_names;
      std::stack<TrigConf::TriggerItemNode*>nodes;
      TrigConf::TriggerItemNode*node;

      threshold_names.clear();
      node = item->topNode();
      nodes.push( item->topNode() );

      // collect unique list (= set) of threshold names for this item
      while (!nodes.empty()) {
         node = nodes.top(); nodes.pop();
         if (node == NULL)
            continue;
         if (node->isThreshold()) {
            if (node->triggerThreshold()) {
               // available if thresholds have been read in
               if (!node->triggerThreshold()->name().empty())
                  threshold_names.insert(node->triggerThreshold()->name());
            } else if (!node->thresholdName().empty()) {
               // fall back solution
               threshold_names.insert(node->thresholdName());
            }
         }
         BOOST_FOREACH(TrigConf::TriggerItemNode* childnode, node->children()) {
            nodes.push(childnode);
         }
      }

      // collect corresponding TEs and add them using appendFeatures()
      tes.clear();
      tes.resize(threshold_names.size());   
      tesit = tes.begin();
      std::set< std::string >::iterator setstrit;

      for (setstrit = threshold_names.begin(); setstrit != threshold_names.end(); ++setstrit, ++tesit) {
         cgm()->navigation()->getAllOfType(TrigConf::HLTUtils::string2hash(*setstrit), *tesit, true);
      }

      appendFeatures(tes, f);    
   }
    
   if (msgLvl(MSG::DEBUG))
      log() << MSG::DEBUG << "features: features container size: "<< f.getCombinations().size() << endreq;
   return f;
}


void Trig::ChainGroup::appendFeatures(std::vector< std::vector< HLT::TriggerElement*> >& tes, Trig::FeatureContainer& fc) const {

  // appends (combinations of) TriggerElements to FeatureContainer
  if (tes.size() == 0) // ComboIterator::isValid would return true in this case
    return;

  HLT::ComboIterator combination(tes, cgm()->navigation());
  while (combination.isValid()) {
    fc.addWithChecking(Combination(*combination, cgm())); 

    if (msgLvl(MSG::VERBOSE))
      log() << MSG::VERBOSE << " adding combination" << Combination(*combination, cgm()) << endreq;

    ++combination;
  }
}


namespace ChainGroup_impl {
  bool vsize(const std::vector<TrigConf::HLTTriggerElement*>& a , const std::vector<TrigConf::HLTTriggerElement*>& b) {
    return a.size() < b.size();
  }
}
/*
bool Trig::ChainGroup::getTEs(const TrigConf::HLTChain* conf, std::vector< std::vector<HLT::TriggerElement*> >& tes,
                              unsigned int condition) const
{
  if ( condition != TrigDefs::Physics && condition != TrigDefs::alsoDeactivateTEs ) {
    throw std::runtime_error("Only two flags can be supplied to features");
  }

  std::vector<TrigConf::HLTSignature*> sigs  = conf->signatureList();

  if (sigs.size() == 0 )
    return true;  // this is chain w/o signatures
  if ( !sigs.back()) 
    return false; // this is bad condition, we are missing configuration
 
  const std::vector<TrigConf::HLTTriggerElement*>& finalTypes = sigs.back()->outputTEs();
  if ( log().level() <= MSG::VERBOSE ) {
    log() << MSG::VERBOSE << "getTEs: finalTEs list: " << finalTypes.size() << endreq;
    
  }
  std::vector<unsigned int> finalTEIds;
  std::back_insert_iterator<std::vector<unsigned int> > idsInserter(finalTEIds);
  transform(finalTypes.begin(), finalTypes.end(), idsInserter, std::mem_fun(&TrigConf::HLTTriggerElement::id));

  if ( log().level() <= MSG::VERBOSE ) {
    log() << MSG::VERBOSE << "getTEs: finalTEs ids: " << finalTEIds << endreq;
  }
  std::vector< std::vector<TrigConf::HLTTriggerElement*> > teTypes;

  // in case we only need passing TEs we call this method with second arg beeing true
  collectConfigurationTEs(conf, condition == TrigDefs::Physics, teTypes);
  log() << MSG::VERBOSE << "getTEs: all configuration TEs: " << teTypes << endreq;

  // we know now all the types
  // they are in 2D array 
  // 1st dimension is number of chain steps
  // 2nd dimension is number of TEs in signatures
  // now we need to loop over number of dimensions of the and collect terminal TEs
  
  unsigned int sigSize  = 0;
  if ( !teTypes.empty() ) {
    // calculate longest signature 
    sigSize = max_element(teTypes.begin(), teTypes.end(), ChainGroup_impl::vsize)->size(); // find max
  }
  if ( log().level() <= MSG::VERBOSE ) {
    log() << MSG::VERBOSE << "getTEs: sig size: " << sigSize << endreq;
  }

  for ( unsigned int teI = 0; teI < sigSize; ++teI ) {
    std::vector<HLT::TriggerElement*> terminalTEs;
    for ( unsigned int sigI = 0; sigI < teTypes.size(); ++sigI ) {
      // protect against shorter signatures
      if ( teI < teTypes[sigI].size() )  
	collectNavigationTEs(teTypes[sigI][teI]->id(), (condition==TrigDefs::Physics), finalTEIds, terminalTEs);
    }
    if ( log().level() <= MSG::VERBOSE ) {
      log() << MSG::VERBOSE << "getTEs: terminalTEs size: " << terminalTEs.size() << endreq;
    }
    tes.push_back(terminalTEs);
  }
  return true;
}


bool Trig::ChainGroup::collectConfigurationTEs(const TrigConf::HLTChain* conf, bool onlyFinals, 
					       std::vector< std::vector<TrigConf::HLTTriggerElement*> >& tes) const {
  std::vector<TrigConf::HLTSignature*> sigs  = conf->signatureList();
  if (sigs.size() == 0 || !sigs.back()) return true;

  if (onlyFinals) {
    tes.push_back( sigs.back()->outputTEs());

  } else {
    //  const std::vector<TrigConf::HLTTriggerElement*>& finalTypes = sigs.back()->outputTEs();      
    std::vector<TrigConf::HLTSignature*>::const_iterator typesIt;
    for ( typesIt = sigs.begin(); typesIt != sigs.end(); ++typesIt ) {
      tes.push_back( (*typesIt)->outputTEs() );
    }
  }
  if ( log().level() <= MSG::VERBOSE ) {
    std::string message;
    std::vector< std::vector<TrigConf::HLTTriggerElement*> >::const_iterator vvIt;
    std::vector<TrigConf::HLTTriggerElement*>::const_iterator vIt;
    for ( vvIt = tes.begin(); vvIt != tes.end(); ++vvIt ) {
      for ( vIt = vvIt->begin(); vIt != vvIt->end(); ++vIt ) {
	message += (*vIt)->name() + " ";
      }
      message += " / "; 
    }
    log() << MSG::VERBOSE << "Config TE: " << message << endreq;
  }

  return true;
}





struct NotLastInThisLevelNorTerminal {
       NotLastInThisLevelNorTerminal(const std::vector<unsigned int>& ids) 
    : m_ids(ids){}
  bool operator()(const HLT::TriggerElement* te) {
    if (  HLT::NavigationCore::isTerminalNode(te) )
      return false;
    if ( count(m_ids.begin(), m_ids.end(), te->getId()))
      return false;
    return true;
  }
  const std::vector<unsigned int>&  m_ids;
};


bool Trig::ChainGroup::collectNavigationTEs(unsigned int TEid,
					    bool activeOnly,
					    const std::vector<unsigned int>& lastIds,
					    std::vector<HLT::TriggerElement*>& tes) const {
  std::vector<HLT::TriggerElement*> thisTypeTEs;
  std::back_insert_iterator<std::vector<HLT::TriggerElement*> > outputIt( tes );

  cgm()->navigation()->getAllOfType(TEid, thisTypeTEs, activeOnly);
  if ( log().level() <= MSG::VERBOSE ) {
    log() << MSG::VERBOSE << "collectNavigationTEs: id: " << TEid
	  << " got back " << thisTypeTEs 
	  << " active " << activeOnly << endreq;
  }
  if ( !activeOnly )
    remove_copy_if(thisTypeTEs.begin(), thisTypeTEs.end(), outputIt, NotLastInThisLevelNorTerminal(lastIds));
  else
    copy(thisTypeTEs.begin(), thisTypeTEs.end(), outputIt);

  return true;
}

*/	
