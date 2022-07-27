/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
#ifndef TRIGSERVICES_TRIGMESSAGESVC_H
#define TRIGSERVICES_TRIGMESSAGESVC_H

// Include files
#include <iosfwd>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include "tbb/concurrent_queue.h"

#include <TH1I.h>
#include <TH2I.h>

#include "CxxUtils/checker_macros.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/Message.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/StatusCode.h"

// Helper to mark some virtual methods as not supported
#define NOTSUPPORTED                                                                               \
  throw std::logic_error(std::string(__func__) + " is not supported by TrigMessageSvc")

// Forward declarations
class ISvcLocator;
class TH1I;
class TH2I;

/**@class TrigMessageSvc
 * @brief MessageSvc used by the HLT applications
 *
 * This MessageSvc implementation it used by the HLT applications. It has some additional
 * features compared to the default Gaudi MessageSvc. Most notably the forwarding of messages
 * to the TDAQ ERS message system.
 *
 * The message suppression is configured with the following parameters:
 * @param <level>Limit = 0:       no message suppression for \<level\>
 * @param <level>Limit = N > 0:   suppress messages after N messages (per source)
 * @param <level>Limit = -N < 0:  use logarithmic suppression after N messages (per message)
 *
 * Note, that the logarithmic suppression works on a per-message basis (ignoring any digits
 * in the message).
 *
 * @author    Iain Last, Werner Wiedenmann, Frank Winklmeier
 */
class TrigMessageSvc : public extends<Service, IMessageSvc, IIncidentListener> {
public:
  typedef std::map<std::string, int, std::less<> > ThresholdMap;

  TrigMessageSvc(const std::string& name, ISvcLocator* svcloc);

  virtual StatusCode reinitialize() override;
  virtual StatusCode initialize() override;
  virtual StatusCode start() override;
  virtual StatusCode stop() override;
  virtual StatusCode finalize() override;
  virtual void handle( const Incident& incident ) override;

  virtual void reportMessage(const Message& message) override;
  virtual void reportMessage(const Message& msg, int outputLevel) override;
  virtual void reportMessage(std::string source, int type, std::string message) override;
  virtual std::ostream* defaultStream ATLAS_NOT_CONST_THREAD_SAFE() const override
  {
    return m_defaultStream;
  }

  virtual void setDefaultStream(std::ostream* stream) override
  {
    m_defaultStream = stream;
  }

  virtual int outputLevel() const override;
  virtual int outputLevel(std::string_view source) const override;
  virtual void setOutputLevel(int new_level) override;
  virtual void setOutputLevel(std::string_view source, int new_level) override;
  virtual int messageCount(MSG::Level logLevel) const override;

  virtual bool useColor() const override { return m_color; }
  virtual std::string getLogColor(int) const override { return ""; }

  ///@{ Not supported by this implementation
  virtual void reportMessage(const StatusCode&, std::string_view) override { NOTSUPPORTED; }
  virtual void insertMessage(const StatusCode&, Message) override { NOTSUPPORTED; }
  virtual void eraseMessage() override { NOTSUPPORTED; }
  virtual void eraseMessage(const StatusCode&) override { NOTSUPPORTED; }
  virtual void eraseMessage(const StatusCode&, const Message&) override { NOTSUPPORTED; }
  virtual void insertStream(int, std::string, std::ostream*) override { NOTSUPPORTED; }
  virtual void eraseStream() override { NOTSUPPORTED; }
  virtual void eraseStream(int) override { NOTSUPPORTED; }
  virtual void eraseStream(int, std::ostream*) override { NOTSUPPORTED; }
  virtual void eraseStream(std::ostream*) override { NOTSUPPORTED; }
  ///@}

private:
  //////////////////////////////////////////////////////////////////////
  // Properties
  //////////////////////////////////////////////////////////////////////
  Gaudi::Property<std::string> m_defaultFormat{this, "Format", Message::getDefaultFormat(),
                                               "Default message format"};
  Gaudi::Property<std::string> m_ersFormat{this, "ErsFormat", Message::getDefaultFormat(),
                                           "ERS message format"};
  Gaudi::Property<std::string> m_defaultTimeFormat{
      this, "timeFormat", Message::getDefaultTimeFormat(), "Message time format"};
  Gaudi::Property<bool> m_stats{this, "showStats", false, "Show message statistics"};
  Gaudi::Property<unsigned int> m_statLevel{this, "statLevel", 0,
                                            "Show total message statistics for >= level"};
  Gaudi::Property<unsigned int> m_publishLevel{this, "publishLevel", MSG::INFO,
                                               "Publish message statistics for this and higher message levels"};
  Gaudi::Property<unsigned int> m_eventIDLevel{this, "printEventIDLevel", MSG::NIL,
                                               "Print event ID for this and higher message levels"};
  Gaudi::Property<bool> m_color{this, "useColors", false,
                                "Colors are not supported by TrigMessageSvc"};
  Gaudi::Property<bool> m_suppress{this, "enableSuppression", false, "Enable message suppression"};
  Gaudi::Property<bool> m_suppressRunningOnly{this, "suppressRunningOnly", true,
                                              "Use message suppression only during RUNNING state"};

  std::array<Gaudi::Property<std::vector<std::string>>, MSG::NUM_LEVELS> m_thresholdProp{
      {{/*ignored*/},
       {this, "setVerbose"},
       {this, "setDebug"},
       {this, "setInfo"},
       {this, "setWarning"},
       {this, "setError"},
       {this, "setFatal"},
       {this, "setAlways"}}};

  std::array<Gaudi::Property<int>, MSG::NUM_LEVELS> m_msgLimit{{{this, "defaultLimit", 500},
                                                                {this, "verboseLimit", 500},
                                                                {this, "debugLimit", 500},
                                                                {this, "infoLimit", 500},
                                                                {this, "warningLimit", 500},
                                                                {this, "errorLimit", 500},
                                                                {this, "fatalLimit", 500},
                                                                {this, "alwaysLimit", 0}}};

  /**
   * Special properties to control output to ERS of individual sources.
   * The syntax is as follows (these are NOT regular expressions):
   *
   * useErsFatal = []                       # forward none (default)
   * useErsFatal = ['*']                    # forward all
   * useErsFatal = ['CoreDumpSvc','MyAlg']  # forward these sources
   * useErsFatal = ['*','!MyAlg']           # forward all except MyAlg
   */
  std::array<Gaudi::Property<std::vector<std::string>>, MSG::NUM_LEVELS> m_useERS{
      {{/*ignored*/},
       {this, "useErsVerbose", {}},
       {this, "useErsDebug", {}},
       {this, "useErsInfo", {}},
       {this, "useErsWarning", {}},
       {this, "useErsError", {}},
       {this, "useErsFatal", {}},
       {this, "useErsAlways", {}}}};

  Gaudi::Property<int> m_ersEventLimit{this, "ersPerEventLimit", -1,
     "Maximum number of messages (per event and level) that are forwarded to ERS (-1: disabled)"};

  //////////////////////////////////////////////////////////////////////
  // Private members
  //////////////////////////////////////////////////////////////////////
  std::ostream* m_defaultStream = &std::cout; ///< Pointer to the output stream.
  ThresholdMap m_thresholdMap;                ///< Output level threshold map

  /// Private helper class to keep the count of messages of a type (MSG::LEVEL).
  struct MsgAry final {
    /// Internal array of counters.
    std::array<int, MSG::NUM_LEVELS> msg = {{0}};
    /// Default constructor.
    MsgAry() = default;
  };

  std::map<std::string, MsgAry> m_sourceMap;     ///< counts per source
  std::array<int, MSG::NUM_LEVELS> m_msgCount{};   ///< counts per level
  std::map<size_t, unsigned int> m_msgHashCount; ///< counts per message hash
  std::unordered_map<EventContext::ContextID_t,
                     std::pair<EventContext::ContextEvt_t, MsgAry>> m_slotMsgCount; ///< counts per slot and level

  bool m_doPublish{false};  ///< are we publishing message statistics?
  bool m_doSuppress{false}; ///< is suppression currently enabled?

  TH1I* m_msgCountHist{nullptr};    ///< Message counting per level histogram
  TH2I* m_msgCountSrcHist{nullptr}; ///< Message counting per message source

  mutable std::recursive_mutex m_thresholdMapMutex; /// (@see MsgStream::doOutput).

  bool m_asyncReporting{false};     ///< Async reporting active
  std::thread m_thread;             ///< Thread for asynchronous reporting
  tbb::concurrent_bounded_queue<std::function<void()>> m_messageActionsQueue;

  void setupLimits(Gaudi::Details::PropertyBase& prop);
  void setupThreshold(Gaudi::Details::PropertyBase& prop);
  bool passErsFilter(const std::string& source, const std::vector<std::string>& filter) const;
  bool passErsLimit(const Message& msg);
  void i_reportMessage(const Message& msg, int outputLevel);
  void i_reportERS(const Message& msg) const;
  void asyncReporting();
  void bookHistograms();
};

#endif
