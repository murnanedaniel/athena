// -*- c++ -*-
/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef RNGCOMPS_ATHRNGSVC_H
#define RNGCOMPS_ATHRNGSVC_H

#include <unordered_map>
#include <vector>
#include <string>
#include <mutex>
#include <functional>

#include "AthenaKernel/IAthRNGSvc.h"
#include "AthenaBaseComps/AthService.h"
#include "GaudiKernel/IIncidentListener.h"

namespace ATHRNG {
  class RNGWrapper;
}
namespace CLHEP{
  class HepRandomEngine;
}

/// @class AthRNGSvc
/// @brief A service to manage multiple RNG streams in thread-safe way.
///
/// The random engines are provided via the RNGWrapper which dereferences to
/// the appropriate slot-local engine.
///
/// @todo Move from manual pointer management to smart pointer management.
///
class AthRNGSvc : public extends<AthService, IAthRNGSvc, IIncidentListener>
{

public:
  virtual ATHRNG::RNGWrapper* GetEngine(const std::string& streamName) override final;
  virtual ATHRNG::RNGWrapper* setOnDefinedSeeds(uint64_t theSeed,
                                                const std::string& streamName) override final;
  virtual ATHRNG::RNGWrapper* setOnDefinedSeeds(uint64_t eventNumber,
                                                uint64_t runNumber,
                                                const std::string& streamName) override final;
  ///seed all streams we manage, combining theSeed and the stream names
  virtual bool setAllOnDefinedSeeds (uint64_t theSeed) override final;
  ///seed all streams, combining eventNumber, runNumber and the stream names
  virtual bool setAllOnDefinedSeeds (uint64_t eventNumber, uint64_t runNumber) override final;
  AthRNGSvc(const std::string& name, ISvcLocator* svc);
  virtual ~AthRNGSvc();
  virtual void handle( const Incident& incident );
  virtual void print(const std::string& streamName) override final;
  virtual void print() override final;
  StatusCode initialize() override final;
  StatusCode start() override final;
  StatusCode finalize() override final;

private:
  void CreateStream(uint64_t seed1, uint64_t seed2,
                    const std::string& streamName);
  size_t hashCombine(size_t h1, size_t h2){ return (h1^(h2+(h1<<6)+(h1>>2)));};
  std::unordered_map<std::string,
    std::pair<ATHRNG::RNGWrapper*, CLHEP::HepRandomEngine*>> m_wrappers;
  std::mutex m_mutex;
  std::string m_RNGType;
  std::size_t m_numSlots;
  std::vector<std::string> m_seeds;
  bool m_initialized;
  typedef std::function<CLHEP::HepRandomEngine*(void)> factoryFunc;
  factoryFunc m_fact;
};

#endif
