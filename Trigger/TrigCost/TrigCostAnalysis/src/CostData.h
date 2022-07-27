/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRIGCOSTANALYSIS_COSTDATA_H
#define TRIGCOSTANALYSIS_COSTDATA_H 1

#include "GaudiKernel/StatusCode.h"
#include "xAODTrigger/TrigCompositeContainer.h"
#include "TrigConfData/HLTChain.h"
#include "TrigCompositeUtils/AlgToChainTool.h"

#include <map>
#include <vector>

/**
 * @class CostData
 * @brief Caches and propagates event data to be used by monitoring algorithms.
 *
 * Cost Monitors and their Counters need access to different storegate collections, derived data,
 * and other normalisation factors. The CostData object wraps all of these, providing a single
 * "data" object which is passed down to all clients.
 */
class CostData {
  public:
    /**
     * @brief Construct an empty CostData.
     */
    CostData();

    /**
     * @brief Default destructor.
     */
    ~CostData() = default;

    /**
     * @brief Forbid assignment.
     */
    CostData& operator=(const CostData&) = delete;

    /**
     * @brief Forbid copy.
     */
    CostData(const CostData&) = delete;

    /**
     * @brief Cache the cost and ros collections, after formally requesting it from storegate.
     */
    StatusCode set(const xAOD::TrigCompositeContainer* costCollection, const xAOD::TrigCompositeContainer* rosCollection, uint32_t onlineSlot);

    /**
     * @brief Getter of the cached algorithm cost collection pointer.
     */
    const xAOD::TrigCompositeContainer& costCollection() const;

    /**
     * @brief Getter of the cached ros cost collection pointer.
     */
    const xAOD::TrigCompositeContainer& rosCollection() const;

    /**
     * @brief Getter of the ROS to ROB map.
     */
    const std::map<std::string, std::vector<uint32_t>>& rosToRobMap() const;

    /**
     * @brief Set ROS to ROB map
     */
    void setRosToRobMap(const std::map<std::string, std::vector<uint32_t>>& rosToRobMap);

    /**
     * @brief Getter of the alg name to chains map.
     */
    const std::map<std::string, std::set<size_t>>& chainToAlgMap() const;

    /**
     * @brief Set the alg name to chains map.
     */
    void setChainToAlgMap( const std::map<std::string, std::set<size_t>>& algToChains );

    /**
     * @brief Getter of the chain to its unique alg names map.
     */
    const std::map<std::string, std::set<size_t>>& chainToUniqAlgMap() const;

    /**
     * @brief Set the chain to its unique alg names map.
     */
    void setChainToUniqAlgMap( const std::map<std::string, std::set<size_t>>& algToChains );

    /**
     * @brief Getter of the sequence to alg idx map.
     */
    const std::map<std::string, std::map<int16_t, std::set<size_t>>>& sequencersMap() const;

    /**
     * @brief Set the sequence to alg idx map.
     */
    void setSequencersMap( const std::map<std::string, std::map<int16_t, std::set<size_t>>>& seqToAlg );

    /**
     * @brief Getter of the seeded chains set.
     */
    const std::vector<TrigCompositeUtils::AlgToChainTool::ChainInfo>& seededChains() const;

    /**
     * @brief Set the seeded chains set.
     */
    void setSeededChains(const std::vector<TrigCompositeUtils::AlgToChainTool::ChainInfo>& seededChains);

    /**
     * @brief Getter of map between algorithm (index in costCollection) and ROS requests (indicies in rosCollection)
     */
    const std::map<size_t, std::vector<size_t>>& algToRequestMap() const;

    /**
     * @brief Setter of effective P1 walltime represented by the current event.
     */
    void setLb(uint32_t lb);

    /**
     * @brief Setter of the online Slot number of the current event.
     */
    void setOnlineSlot(uint32_t slot);

    /**
     * @brief Setter of effective P1 walltime represented by the current event, or the current lumi block. As specified by the second parameter
     */
    void setLivetime(float time, bool liveTimeIsPerEvent);

    /**
     * @brief Getter of effective P1 walltime represented by either the current event, or the current lumi block.
     * @return Walltime in seconds.
     */
    float liveTime() const;

    /**
     * @brief If a call to liveTime() is providing data on a single event or a whole LB
     * @return Walltime in seconds.
     */
    bool liveTimeIsPerEvent() const;

    /**
     * @brief Current luminosity block number
     * @return Luminosity block number.
     */
    uint32_t lb() const;

    /**
     * @return Online slot number
     */
    uint32_t onlineSlot() const;

    /**
     * @return True if event was processed in the master slot (0), and hence contains cost data spanning all concurrent slots
     */
    bool isMasterSlot() const;

    /**
     * @brief Getter of the total algorithm CPU time in the event. 
     * @return Total CPU time in milliseconds.
     */
    float algTotalTimeMilliSec() const;


    /**
     * @brief Get the class typename given an algorithm instance name. Name is supplied in serialised hashed form.
     */
    const std::string& algNameToClassType(size_t algNameHash) const;

    /**
     * @brief Set internal type map pointer
     * @return Total CPU time in milliseconds.
     */
    void setTypeMap( const std::unordered_map<uint32_t, std::string>& typeMap );

  private:

    /**
     * @brief Compute and cache derived quantities, called automatically after set().
     * Computes algTotalTimeMilliSec()
     */
    StatusCode cache();

    const xAOD::TrigCompositeContainer* m_costCollection; //!< Cached non-owning pointer to main algorithm cost collection.
    const xAOD::TrigCompositeContainer* m_rosCollection = nullptr; //!< Cached non-owning pointer to ros cost collection.
    uint64_t m_algTotalTime; //!< Integrated CPU time of all algorithms in the event. Stored in discrete microseconds.
    float m_liveTime; //!< Effective walltime of either the event or the LB, in seconds (@see m_liveTimeIsPerEvent).
    uint32_t m_lb; //!< Current luminosity block number
    uint32_t m_slot; //!< Current online slot number
    bool m_liveTimeIsPerEvent; //!< If the livetime represents a single event or all of the current LB
    const std::unordered_map<uint32_t, std::string>* m_typeMapPtr; //!< Cached non-owning pointer mapping algorithm instance names to types
    std::map<size_t, std::vector<size_t>> m_algToRos; //!< Mapping of indexes from m_costCollection to corresponding ROS requests made by algorithm
    const std::map<std::string, std::vector<uint32_t>>* m_rosToRob = nullptr; //!< Mapping of ROS corresponding to ROB requests
    const std::map<std::string, std::set<size_t>>* m_chainToAlgIdx = nullptr; //!<Mapping of chain to algorithms idx
    const std::map<std::string, std::set<size_t>>* m_chainToUniqAlgIdx = nullptr; //!<Mapping of chain name to its unique algorithms
    const std::map<std::string, std::map<int16_t, std::set<size_t>>>* m_sequencers = nullptr; //!<Mapping of sequence to algorithms
    const std::vector<TrigCompositeUtils::AlgToChainTool::ChainInfo>* m_seededChains = nullptr; //!<Set of seeded chains to monitor

};

#endif // TRIGCOSTANALYSIS_COSTDATA_H
