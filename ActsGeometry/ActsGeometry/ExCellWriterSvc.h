#ifndef GEOMACTS_EXCELLWRITER_SERVICE
#define GEOMACTS_EXCELLWRITER_SERVICE

#include "AthenaBaseComps/AthService.h"
#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/Property.h"  /*no forward decl: typedef*/

#include "ActsGeometry/IExCellWriterSvc.h"
#include "ActsGeometry/Extrapolation/RootExCellWriter.hpp"

#include "Acts/EventData/TrackParameters.hpp"

#include <vector>
#include <deque>
#include <mutex>
#include <thread>
#include <atomic>


namespace Acts {

template <class>
class ExtrapolationCell;


class ExCellWriterSvc : public extends<AthService, IExCellWriterSvc> {
public:
    
  virtual StatusCode initialize() override;
  virtual StatusCode finalize() override;
    
  ExCellWriterSvc( const std::string& name, ISvcLocator* svc );

  void
  store(std::vector<Acts::ExtrapolationCell<Acts::TrackParameters>>& ecells) override;

private:
  using ExCellCharged = Acts::ExtrapolationCell<Acts::TrackParameters>;

  using queue_item_t = std::pair<size_t, ExCellCharged>;

  std::shared_ptr<RootExCellWriter<Acts::TrackParameters>> m_rootEccWriter;
  std::deque<queue_item_t> m_queue;
  std::mutex m_chargedMutex;
  std::thread m_writeThread;
  std::atomic<bool> m_doEnd;

  void doWrite();

  // jobOptions properties
  Gaudi::Property<std::string> m_filePath{this, "FilePath", "excells_charged.root", "Output root file for charged particle"};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "extrapolation_charged", ""};
  Gaudi::Property<bool> m_writeBoundary{this, "WriteBoundary", true, ""};
  Gaudi::Property<bool> m_writeMaterial{this, "WriteMaterial", true, ""};
  Gaudi::Property<bool> m_writeSensitive{this, "WriteSensitive", true, ""};
  Gaudi::Property<bool> m_writePassive{this, "WritePassive", true, ""};


};

}


#endif 
