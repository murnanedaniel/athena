/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGGER_DECISION_TOOL_CHAIN_GROUP_H
#define TRIGGER_DECISION_TOOL_CHAIN_GROUP_H

/**********************************************************************************
 * @Project: TrigDecisionTool
 * @Package: TrigDecisionTool
 * @class  : ChainGroup
 *
 * @brief container to hold trigger chains
 *
 * @author Michael Begel  <michael.begel@cern.ch>  - Brookhaven National Laboratory
 * @author Joerg Stelzer  <Joerg.Stelzer@cern.ch>  - DESY
 * @author Tomasz Bold     <Tomasz.Bold@cern.ch>   - UC Irvine - AGH-UST Krakow
 *
 ***********************************************************************************/

#include <iterator>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <stack>
#include <boost/algorithm/string.hpp>

#include "TrigDecisionInterface/Conditions.h"
#include "TrigDecisionInterface/GroupProperties.h"
#include "TrigDecisionTool/CacheGlobalMemory.h"
#include "TrigDecisionTool/FeatureContainer.h"
#include "TrigDecisionTool/Logger.h"
#include "TrigDecisionTool/FeatureRequestDescriptor.h"
#include "TrigSteeringEvent/Enums.h"
#include "TrigCompositeUtils/TrigCompositeUtils.h"

namespace HLT {
  class Chain;
}

namespace LVL1CTP {
  class Lvl1Item;
}

namespace TrigConf {
  class TriggerItem;
  class HLTChain;
  class HLTTriggerElement;
}


namespace Trig {

   class ChainGroup : public virtual Logger {

     friend class CacheGlobalMemory;

   public:

      ChainGroup(const std::vector< std::string >& triggerNames,
                 Trig::CacheGlobalMemory& parent);
      ~ChainGroup() = default;

      typedef std::vector<std::string>::const_iterator const_iterator;
      const Trig::ChainGroup& operator+(const Trig::ChainGroup& rhs);
      bool operator==(const Trig::ChainGroup& rhs);
      bool operator!=(const Trig::ChainGroup& rhs);

      /**
       * @brief adds alias (sort understandabel name) to the group
       **/
      void addAlias(const std::string& alias);


      /**
       * @brief tells if chain group passed
       * @param conditions is modifying the question
       **/
      bool isPassed(unsigned int condition=TrigDefs::Physics) const;
      

      /**
       * @brief returns prescale factor
       * for chain group with single chain in returns real prescale factor
       * for real chain group composed of many prescaled chains returns 1 if at least one chain is unprescaled
       **/
      float getPrescale(unsigned int condition=TrigDefs::Physics) const;


      std::vector< std::string > getListOfTriggers() const;
      std::vector< std::string > getListOfStreams() const;
      std::vector< std::string > getListOfGroups() const;
      std::vector< std::string > getListOfThresholds() const;
      std::vector< std::string > getListOfSignatures() const;
      std::vector< std::vector< std::string > > getListOfTriggerElements() const;
      std::vector< std::vector< TrigConf::HLTTriggerElement* > > getHLTTriggerElements() const;

      
      /**
       * @brief returns bits (OR ed) of the chain group
       * Meaning of the returned bits can be understood by using masks defined in TrigDefs
       **/
      unsigned int isPassedBits() const;


      /**
       * @brief returns most severe error in the chains composing that chain group
       * for L1 it is just OK
       * If there is suspicion that there are other problems in the CG one needs to loop over chains and 
       * check each of them.
       **/
      HLT::ErrorCode error() const;


      /**
       * @brief returns all features related to given chain group of HLT chains or L1 items
       * Note: This does not yet work for L1_FJ..., i.e. no features are returned for these items.
       **/
      const FeatureContainer features(unsigned int condition = TrigDefs::Physics) const;
    
      /**
       * @brief returns typed features related to given chain group of HLT chains or L1 items
       * Note: This is a RUN 3 (and on) function.
       * @param[in] eventStore Event store pointer. To migrate to readHandles with the rest of the TDT soon
       * @param[in] HLTSummaryKeyIn SG Key to the navigation summary container
       * @param[in] condition Condition requirement. Only TrigDefs::Physics and TrigDefs::includeFailedDecisions are supported.
       * @param[in] containerSGKey Optional requirement to return only features within the specified container name. Not checked if not specified. 
       * @param[in] featureCollectionMode For lastFeatureOfType, stop exploring each route through the navigation once one matching feature has been found.
       * @param[in] navElementLinkKey Optional name of element link as saved online. The "feature" link is enforced, others may have been added. 
       * @param[in] restrictToLegIndex Optional index of a leg for mult-leg chains. Features will only be returned on the specified leg. Default is all legs.
       * @return Vector of LinkInfo, where each entry wraps an ElementLink to the feature, and the Decision object it came from.
       **/  
      template<class CONTAINER>
      std::vector< TrigCompositeUtils::LinkInfo<CONTAINER> > features(const asg::EventStoreType* eventStore,
                const SG::ReadHandleKey<TrigCompositeUtils::DecisionContainer>& HLTSummaryKeyIn,
                unsigned int condition = TrigDefs::Physics,
                const std::string& containerSGKey = "",
                const unsigned int featureCollectionMode = TrigDefs::lastFeatureOfType,
                const std::string& navElementLinkKey = TrigCompositeUtils::featureString(),
                const int          restrictToLegIndex = -1) const;

      // 
      const std::vector< std::string >& patterns() const {return m_patterns;}
   private:

      bool  isCorrelatedL1items(const std::string& item) const;
      float correlatedL1Prescale(const std::string& item) const;
      float calculatePrescale(unsigned int condition=TrigDefs::Physics);

      void appendFeatures(std::vector< std::vector< HLT::TriggerElement*> >& tes, FeatureContainer& fc) const;

      /**
       * @brief names of triggers within chain group
       **/
      const std::vector< std::string >& names() const {return m_names;}

      bool HLTResult(const std::string& chain, unsigned int condition) const;
      bool L1Result(const std::string& item, unsigned int condition) const;

      unsigned int HLTBits(const std::string& chain, const std::string& level, const TrigCompositeUtils::DecisionIDContainer& passExpress) const;
      unsigned int L1Bits(const std::string& item) const;

      float HLTPrescale(const std::string& chain, unsigned int condition) const;
      float L1Prescale(const std::string& item, unsigned int condition) const;

      std::string getLowerName(const std::string& EFname) const;

      std::vector<std::string> m_patterns;  //!< patterns with which the CG was constructed
    
      std::set<const TrigConf::HLTChain*>           m_confChains;
      std::set<const TrigConf::TriggerItem*>        m_confItems;

#ifndef __REFLEX__
      // quick cache (external therefore reference) of the result per event
      Trig::CacheGlobalMemory&                       m_cgm;
#endif
      std::vector< std::string > m_names; //!< names of trigger derived from patterns & current configuration

      const Trig::CacheGlobalMemory& cgm_assert() const;
      const Trig::CacheGlobalMemory& cgm() const { return m_cgm; }
      Trig::CacheGlobalMemory& cgm() { return m_cgm; }

      // update the configuration
      void update(const TrigConf::HLTChainList* confChains,
                  const TrigConf::ItemContainer* confItems,
                  TrigDefs::Group prop = TrigDefs::Group::Default);

      ChainGroup& operator= (const ChainGroup&);

      float m_prescale{0};

   };
} // End of namespace

#include "ChainGroup.icc"

#endif
