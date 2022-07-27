/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file   CoreDumpSvc.cxx
 * @brief  Implementation of the CoreDumpSvc
 * @author Frank Winklmeier
 *
 */

// System includes
#include <ctime>
#include <cstdio>
#include <fcntl.h>
#include <errno.h>
#include <signal.h>
#include <time.h>
#include <sys/types.h>
#ifndef __APPLE__
#include <sys/sysinfo.h>
#else
#include <mach/task.h> 
#include <mach/mach_init.h>
#include <unistd.h>
#endif

// Package includes
#include "CoreDumpSvc.h"

// ROOT includes
#include "TSystem.h"

// Gaudi includes
#include "Gaudi/Property.h"
#include "GaudiKernel/IAlgorithm.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IAlgContextSvc.h"
#include "GaudiKernel/IAlgExecStateSvc.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/System.h"
#include "GaudiKernel/ConcurrencyFlags.h"
#include "GaudiKernel/EventContext.h"

// Athena includes
#include "AthenaKernel/IAthenaSummarySvc.h"
#include "StoreGate/StoreGateSvc.h"
#include "CxxUtils/SealCommon.h"
#include "CxxUtils/SealSignal.h"
#include "CxxUtils/SealDebug.h"
#include "CxxUtils/read_athena_statm.h"
#include "CxxUtils/checker_macros.h"

namespace {

  const char* const horizLine = "-------------------------------------------------------------------------------------\n";

   void ExitOnInt( int sig, siginfo_t*, void* ) {
      if ( sig == SIGINT ) {
      // called on user ^C
	 std::cout << std::endl;
         std::cerr << "Athena           CRITICAL stopped by user interrupt\n";
	 raise(SIGKILL);
      }
   }

} // unnamed namespace


/**
 * @brief  Signal handler for CoreDumpSvc
 *
 * All information accessible from the signal handler is in this namespace.
 * It carries a pointer to the CoreDumpSvc instance. Therefore no static
 * methods are needed in the CoreDumpSvc to provide a function pointer.
 */
namespace CoreDumpSvcHandler
{
  typedef std::map<int, struct sigaction> SigHandler_t;
  
  SigHandler_t oldSigHandler;         ///< old signal handlers
  bool callOldHandler(true);          ///< forward calls to old handlers?
  bool dumpCoreFile(false);           ///< dump core file on exit?
  bool stackTrace(false);             ///< produce stack trace?
  bool fastStackTrace(false);         ///< produce fast stack trace using CxxUtils/Seal
  CoreDumpSvc* coreDumpSvc(nullptr);  ///< pointer to CoreDumpSvc
  std::ostream* ostr(&std::cout);     ///< stream for printing

  std::ostream& log ATLAS_NOT_THREAD_SAFE () { return *ostr; }  ///< convenience method for logging

  /**
   * Signal handler for the CoreDumpSvc
   */
  void action ATLAS_NOT_THREAD_SAFE ( int sig, siginfo_t *info, void* extra )
  {
    // Careful: don't do anything here that might allocate memory.

    // Protect against recursion.
    // We originally used a thread_local here --- but accessing
    // a thread_local can result in a call to malloc.

    const int maxcalls = 64;
    static std::atomic<int> ncalls (0);
    if (++ncalls >= maxcalls) _exit (98);

    static std::mutex tidlist_mutex;
    static size_t ntids ATLAS_THREAD_SAFE = 0;
    static pthread_t tids[maxcalls] ATLAS_THREAD_SAFE;
    {
      pthread_t self = pthread_self();
      std::lock_guard<std::mutex> lock (tidlist_mutex);
      for (size_t i = 0; i < ntids; i++) {
        if (pthread_equal (self, tids[i])) return;
      }
      if (ntids == maxcalls) _exit (98);
      tids[ntids++] = self;
    }

    // Count the number of threads trying to dump.
    static std::atomic<int> inThreads = 0;
    ++inThreads;

    const unsigned int timeoutSeconds = static_cast<unsigned int>(round(coreDumpSvc->m_timeout * 1e-9));

    if ( sig == SIGALRM) {
      if (dumpCoreFile) {
        log() << "Received SIGALRM. Aborting job..." << std::endl;
        // Restore default abort handler that should create a core file
        Athena::Signal::revert (SIGABRT);
        std::abort();
      }
      else {
        log() << "Received SIGALRM. Terminating job..." << std::endl;
        _exit(97);   // exit without raising any further signals
      }
    }

    // Only allow one thread past at a time.
    // Try to assume as little as possible about the state of the library.
    // We don't want to hang forever here, but we also don't want
    // to call any library functions that might use signals under the hood.
    // So use nanosleep() to do the delay --- that's defined to be
    // independent of signals.
    static std::mutex threadMutex;
    const timespec one_second { 1, 0 };
    {
      unsigned int waits = 0;
      while (!threadMutex.try_lock()) {
        nanosleep (&one_second, nullptr);
        if (++waits > timeoutSeconds) _exit (97);
      }
    }

    // setup timeout
    if ( timeoutSeconds > 0 && (sig == SIGSEGV || sig == SIGBUS || sig == SIGABRT) ) {
      // This will trigger SIGALRM, which we then handle ourselves above
      alarm(timeoutSeconds);
    }

    // Do fast stack trace before anything that might touch the heap.
    // For extra paranoia, avoid iostreams/stdio and use write() directly.
    if (fastStackTrace) {
      write (1, horizLine, strlen(horizLine));
      const char* msg = "Producing (fast) stack trace...\n";
      write (1, msg, strlen (msg));
      write (1, horizLine, strlen(horizLine));
      Athena::Signal::fatalDump (sig, info, extra,
                                 Athena::DebugAids::stacktraceFd(),
                                 Athena::Signal::FATAL_DUMP_SIG +
                                 Athena::Signal::FATAL_DUMP_CONTEXT +
                                 Athena::Signal::FATAL_DUMP_STACK);
      write (1, "\n", 1);
    }

    std::cout.flush();
    std::cerr.flush();
    
    if (coreDumpSvc) {
      coreDumpSvc->setSigInfo(info);
      coreDumpSvc->print();
    }

    if (gSystem && stackTrace) {
      log() << horizLine << "Producing stack trace (can be slow, check gdb process)...\n"
            << horizLine << std::flush;
      gSystem->StackTrace();
      log() << std::endl;
    }

    if (callOldHandler) {
      // Call previous signal handler
      // Need to distinguish between the two different types
      const struct sigaction& oact = oldSigHandler[sig];
      log() << horizLine << "Invoking previous signal handler (can be slow, check gdb process)...\n"
            << horizLine << std::flush;
      if ( oact.sa_flags & SA_SIGINFO ) {
        oact.sa_sigaction(sig, info, extra);
      }
      else if (oact.sa_handler != SIG_DFL && oact.sa_handler != SIG_IGN ) {
        oact.sa_handler(sig);
      }
      else {
        log() << "Could not invoke previous signal handler" << std::endl;
      }
    }

    // This thread is done dumping.
    threadMutex.unlock();
    --inThreads;

    if (coreDumpSvc && (sig == SIGSEGV || sig == SIGBUS || sig == SIGABRT) ) {
      // Don't terminate the program while there are other threads
      // trying to dump (but don't wait forever either).
      unsigned int waits = 0;
      while (inThreads > 0 && waits < timeoutSeconds) {
        nanosleep (&one_second, nullptr);
      }

      if (dumpCoreFile) {
        log() << "Aborting job... " << std::endl;
        // Restore default abort handler that should create a core file
        Athena::Signal::revert (SIGABRT);
        std::abort();
      }

      // Exit now on a fatal signal; otherwise, we can hang.
      _exit (99);
    }
  }

}

//================================================================================
// C'tor, D'tor, Property handler
//================================================================================
CoreDumpSvc::CoreDumpSvc( const std::string& name, ISvcLocator* pSvcLocator ) :
  base_class( name, pSvcLocator )
{
  // Set us as the current instance
  CoreDumpSvcHandler::coreDumpSvc = this;
  
  m_callOldHandler.declareUpdateHandler(&CoreDumpSvc::propertyHandler, this);
  m_dumpCoreFile.declareUpdateHandler(&CoreDumpSvc::propertyHandler, this);
  m_stackTrace.declareUpdateHandler(&CoreDumpSvc::propertyHandler, this);
  m_fastStackTrace.declareUpdateHandler(&CoreDumpSvc::propertyHandler, this);
  m_coreDumpStream.declareUpdateHandler(&CoreDumpSvc::propertyHandler, this);
  m_fatalHandlerFlags.declareUpdateHandler(&CoreDumpSvc::propertyHandler, this);
  m_killOnSigInt.declareUpdateHandler(&CoreDumpSvc::propertyHandler, this);
  // Allocate for 2 slots just for now.
  m_usrCoreDumps.resize(2);
  m_sysCoreDumps.resize(2); 
}

CoreDumpSvc::~CoreDumpSvc()
{
  CoreDumpSvcHandler::coreDumpSvc = nullptr;
}

void CoreDumpSvc::propertyHandler(Gaudi::Details::PropertyBase& p)
{
  CoreDumpSvcHandler::callOldHandler = m_callOldHandler;
  CoreDumpSvcHandler::dumpCoreFile = m_dumpCoreFile;
  CoreDumpSvcHandler::stackTrace = m_stackTrace;
  CoreDumpSvcHandler::fastStackTrace = m_fastStackTrace;

  if ( p.name()==m_coreDumpStream.name() ) {
    const std::string val = p.toString();
    if ( val=="stdout" ) {
      CoreDumpSvcHandler::ostr = &std::cout;
    }
    else if ( val=="stderr" ) {
      CoreDumpSvcHandler::ostr = &std::cerr;
    }
    else {
      ATH_MSG_WARNING("'" << val << "' not valid for " << m_coreDumpStream.name()
                      << ": " << m_coreDumpStream.documentation());
    }
  } else if ( p.name() == m_fatalHandlerFlags.name() ) {
    if (m_fatalHandlerFlags.fromString(p.toString()).isSuccess()) {
      if (m_fatalHandlerFlags != 0) {
	Athena::Signal::handleFatal(nullptr, IOFD_INVALID, nullptr, nullptr, m_fatalHandlerFlags);
      }
    } else {
      ATH_MSG_INFO("could not convert [" << p.toString() << "] to integer");
    }
  }
  else if (p.name() ==  m_killOnSigInt.name()) {
    if (m_killOnSigInt.fromString(p.toString()).isSuccess()) {
      if (m_killOnSigInt) {
	ATH_MSG_DEBUG("Will kill job on SIGINT (Ctrl-C)");
	Athena::Signal::handle( SIGINT, ExitOnInt );
      }
    }
    else {
      ATH_MSG_WARNING("Could not convert [" << p.toString() << "] to bool");
    }
  }

}

//================================================================================
// IService implementation
//================================================================================
StatusCode CoreDumpSvc::initialize()
{
  if (m_fatalHandlerFlags != 0) {
      ATH_MSG_INFO("install f-a-t-a-l handler... (flag = " << m_fatalHandlerFlags.value() << ")");
      Athena::Signal::handleFatal(nullptr, IOFD_INVALID, nullptr, nullptr, m_fatalHandlerFlags);
  }

  if (m_killOnSigInt) {
      ATH_MSG_DEBUG("Will kill job on SIGINT (Ctrl-C)");
      Athena::Signal::handle( SIGINT, ExitOnInt );
  }
  
  if ( installSignalHandler().isFailure() ) {
    ATH_MSG_ERROR ("Could not install signal handlers");
    return StatusCode::FAILURE;
  }

  // Register incident handler
  ServiceHandle<IIncidentSvc> incSvc("IncidentSvc", name());
  if ( !incSvc.retrieve().isSuccess() ) {
    ATH_MSG_WARNING ("Unable to retrieve the IncidentSvc");
  }
  else {
    incSvc->addListener(this, IncidentType::BeginRun);
    incSvc->addListener(this, IncidentType::BeginEvent);
    incSvc->addListener(this, IncidentType::EndRun);
    incSvc->addListener(this, IncidentType::EndEvent);
    incSvc->addListener(this,"StoreCleared");
  }
  
  return StatusCode::SUCCESS;
}

StatusCode CoreDumpSvc::start()
{
  auto numSlots = std::max<size_t>(1, Gaudi::Concurrency::ConcurrencyFlags::numConcurrentEvents());
  m_usrCoreDumps.resize(numSlots);
  m_sysCoreDumps.resize(numSlots);
  return StatusCode::SUCCESS;
}

StatusCode CoreDumpSvc::finalize()
{
  ATH_MSG_DEBUG ("Finalizing " << name());

  if ( uninstallSignalHandler().isFailure() ) {
    ATH_MSG_WARNING ("Could not uninstall signal handlers");
    return StatusCode::FAILURE;
  }
   
  return StatusCode::SUCCESS;
}

//================================================================================
// ICoreDumpSvc implementation
//================================================================================

//----------------------------------------------------------------------
// Set a name/value pair in the core dump record
//----------------------------------------------------------------------
void CoreDumpSvc::setCoreDumpInfo( const std::string& name, const std::string& value )
{
  setCoreDumpInfo(Gaudi::Hive::currentContext(), name, value);
}

void CoreDumpSvc::setCoreDumpInfo( const EventContext& ctx, const std::string& name, const std::string& value )
{
  auto slot = ctx.valid() ? ctx.slot() : 0;  
  m_usrCoreDumps.at(slot)[name] = value;
}

//----------------------------------------------------------------------
// Print all core dump records
//----------------------------------------------------------------------
void CoreDumpSvc::print ATLAS_NOT_THREAD_SAFE ()
{
  // Print a FATAL message but don't use the MsgStream anymore once we crashed
  CoreDumpSvcHandler::log() << name() << "   FATAL Caught fatal signal. Printing details to "
                            << m_coreDumpStream.value()
                            << (m_dumpCoreFile ? ". Will try to produce a core dump file on exit." : ".")
                            << std::endl;

  CoreDumpSvcHandler::log() << dump() << std::flush;
}

//----------------------------------------------------------------------
// Print all core dump records
//----------------------------------------------------------------------
std::string CoreDumpSvc::dump() const
{
  std::ostringstream os;
  char buf[26];
  const time_t now = time(nullptr);
  
  os << "-------------------------------------------------------------------------------------" << "\n";
  os << "Core dump from " << name() << " on " << System::hostName()
     << " at " << ctime_r(&now, buf) /*<< "\n"*/; // ctime adds "\n"
  os << "\n";

  // Print additional information if available
  if (m_siginfo) {
    int signo = m_siginfo->si_signo;  // shorthand
    
    os << "Caught signal " << signo
       << "(" << strsignal(signo) << "). Details: "
       << "\n";   

    os << "  errno = " << m_siginfo->si_errno
       << ", code = " << m_siginfo->si_code
       << " (" << Athena::Signal::describe(signo, m_siginfo->si_code) << ")"
       << "\n";
    
    os << "  pid   = " << m_siginfo->si_pid
       << ", uid = " << m_siginfo->si_uid
       << "\n";
    
#ifndef __APPLE__
    // These are set if the POSIX signal sender passed them.
    os << "  value = (" << m_siginfo->si_int << ", "
       << std::hex << m_siginfo->si_ptr << ")" << std::dec << "\n";
#endif

    // memory usage informations
    athena_statm s = read_athena_statm();
    
    const long pagesz = sysconf(_SC_PAGESIZE);
    os << "  vmem = " << s.vm_pages*pagesz/1024./1024.  << " MB\n"
       << "  rss  = " << s.rss_pages*pagesz/1024./1024. << " MB\n";

#ifndef __APPLE__
    // more memory usage informations (system wide stuff)
    // see sysinfo(2)

    {
      struct sysinfo sys;
      if ( 0 == sysinfo(&sys) ) {
        // all sizes are reported in sys.mem_unit bytes
        const float mem_units = sys.mem_unit/(1024.*1024.);
        os << "  total-ram = " << sys.totalram * mem_units << " MB\n"
           << "  free-ram  = " << sys.freeram  * mem_units << " MB\n"
           << "  buffer-ram= " << sys.bufferram* mem_units << " MB\n"
           << "  total-swap= " << sys.totalswap* mem_units << " MB\n"
           << "  free-swap = " << sys.freeswap * mem_units << " MB\n";
      }
    }
#endif

    // This is the interesting address for memory faults.
    if (signo == SIGILL || signo == SIGFPE || signo == SIGSEGV || signo == SIGBUS)
      os << "  addr  = " << std::hex << m_siginfo->si_addr << std::dec << "\n";
    
    os << "\n";
  }
  
  os << "Event counter: " << m_eventCounter << "\n";  


  IAlgExecStateSvc* algExecStateSvc(nullptr);
  IAlgContextSvc* algContextSvc(nullptr);

  // Use AlgExecStateSvc in MT, otherwise AlgContextSvc
  if (Gaudi::Concurrency::ConcurrencyFlags::numConcurrentEvents() > 0) {
    service("AlgExecStateSvc", algExecStateSvc, /*createIf=*/ false).ignore();
  }
  else {
    service("AlgContextSvc", algContextSvc, /*createIf=*/ false).ignore();
  }

  // Loop over all slots
  for (size_t t=0; t < m_sysCoreDumps.size(); ++t){

    // Currently executing algorithm(s)
    std::string currentAlg;
    if (algExecStateSvc) {
      ATH_MSG_DEBUG("Using AlgExecStateSvc to determine current algorithm(s)");
      try {
        // We copy on purpose to avoid modification while we examine it
        auto states = algExecStateSvc->algExecStates(EventContext(0,t));
        for (const auto& kv : states) {
          if (kv.second.state()==AlgExecState::State::Executing)
            currentAlg += (kv.first + " ");
        }
      }
      catch (const GaudiException&) {  // can happen if we get called before any algo execution
        ATH_MSG_INFO("No information from AlgExecStateSvc because no algorithm was executed yet.");
      }
    }
    else if (algContextSvc) {
      ATH_MSG_DEBUG("Using AlgContextSvc to determine current algorithm");
      IAlgorithm* alg = algContextSvc->currentAlg();
      if (alg) currentAlg = alg->name();
    }
    else {
      ATH_MSG_WARNING("AlgExecStateSvc or AlgContextSvc not available. Cannot determine current algorithm.");
    }

    if (currentAlg.empty()) currentAlg = "<NONE>";
    os << "Slot " << std::setw(3) << t << " : Current algorithm = " << currentAlg << std::endl;
        
    // System core dump
    auto &sys = m_sysCoreDumps.at(t);
    if (!sys.LastInc.empty()) {
      os << "         : Last Incident = " << sys.LastInc << std::endl
         << "         : Event ID      = " << sys.EvId << std::endl;
    }
    
    // User core dump
    auto &usr = m_usrCoreDumps.at(t);
    if (!usr.empty()) {
      for (auto &s : usr) {
        os << "         : (usr) " << s.first << " = " << s.second << std::endl;
      }
    }
  }

  if (algContextSvc) {
    os << "Algorithm stack: ";
    if ( algContextSvc->algorithms().empty() ) os << "<EMPTY>" << "\n";
    else {
      os << "\n";
      for (auto alg : algContextSvc->algorithms()) {
        if (alg) os << "   " << alg->name() << "\n";
      }
    }
  }

  os << horizLine;
  os << "| AtlasBaseDir : " << std::setw(66) << getenv("AtlasBaseDir")  << " |\n";
  os << "| AtlasVersion : " << std::setw(66) << getenv("AtlasVersion")  << " |\n";
  os << "| BINARY_TAG   : " << std::setw(66) << getenv("BINARY_TAG")    << " |\n";
  os << horizLine;
  os << " Note: to see line numbers in below stacktrace you might consider running following :\n";
  os << "  atlasAddress2Line --file <logfile>\n";

  IAthenaSummarySvc *iass(nullptr);
  if (service("AthenaSummarySvc",iass,false).isSuccess() && iass) {
    iass->addSummary("CoreDumpSvc",os.str());
    iass->setStatus(1);
    iass->createSummary().ignore();
  }
  
  return os.str();
}

//================================================================================
// IIncidentHandler implementation
//================================================================================

void CoreDumpSvc::handle(const Incident& incident)
{
  //handle is single threaded in context;
  auto slot = incident.context().valid() ? incident.context().slot() : 0;  
  auto &currRec = m_sysCoreDumps.at(slot);

  currRec.LastInc = incident.source() + ":" + incident.type();

  std::ostringstream oss;
  oss << incident.context().eventID();
  currRec.EvId = oss.str();

  if (incident.type()==IncidentType::BeginEvent) {
    // Set up an alternate stack for this thread, if not already done.
    setAltStack();
    ++m_eventCounter;
  } else if (incident.type() == "StoreCleared") {
    // Try to force reallocation.
    auto newstr = currRec.EvId;
    // Intentional:
    // cppcheck-suppress selfAssignment
    newstr[0] = newstr[0];
    currRec.EvId = newstr;
  }

}

//================================================================================
// Helpers for signal handler
//================================================================================

//----------------------------------------------------------------------
// Install signal handler
//----------------------------------------------------------------------
StatusCode CoreDumpSvc::installSignalHandler ATLAS_NOT_THREAD_SAFE ()
{
  ATH_MSG_DEBUG ("Installing signal handler");
  std::ostringstream oss;

  for (auto sig : m_signals) {
#ifndef __APPLE__
    if (sig<1 || sig>SIGRTMAX) {
      ATH_MSG_WARNING ("Invalid signal number " << sig << ". Ignoring.");
      continue;
    }
#endif
    oss << sig << "(" << strsignal(sig) << ") ";

    // Set up an alternate stack for this thread.
    setAltStack();
    
    // Install new signal handler and backup old one
    struct sigaction sigact;
    memset (&sigact, 0, sizeof(sigact));
    sigact.sa_sigaction = CoreDumpSvcHandler::action;
    sigemptyset(&sigact.sa_mask);
    sigact.sa_flags = SA_SIGINFO + SA_ONSTACK;
    int ret = sigaction(sig, &sigact, &(CoreDumpSvcHandler::oldSigHandler[sig]));
    if ( ret!=0 ) {
      ATH_MSG_ERROR ("Error on installing handler for signal " << sig
                     << ": " << strerror(errno));
      return StatusCode::FAILURE;
    }
  }
  ATH_MSG_INFO ("Handling signals: " << oss.str());
  
  return StatusCode::SUCCESS;
}

//----------------------------------------------------------------------
// Uninstall signal handler
//----------------------------------------------------------------------
StatusCode CoreDumpSvc::uninstallSignalHandler ATLAS_NOT_THREAD_SAFE ()
{
  ATH_MSG_DEBUG ("Uninstalling signal handler");

  StatusCode sc = StatusCode::SUCCESS;

  for (const auto& kv : CoreDumpSvcHandler::oldSigHandler) {
    int ret = sigaction(kv.first, &(kv.second), nullptr);
    if ( ret!=0 ) {
      sc = StatusCode::FAILURE;
      ATH_MSG_WARNING("Error on uninstalling handler for signal " << kv.first
                      << ": " << strerror(errno));
    }
  }
  return sc;
}


// Set an alternate stack to use for doing stack traces, so that we
// can continue even if our primary stack is corrupt / exhausted.
// Reserve 2MB on top of the minimum required for a signal handler.
// This sets the alternate stack for the current thread, if it hasn't
// already been done.
void CoreDumpSvc::setAltStack()
{
  std::vector<uint8_t>& stack = s_stack;
  if (stack.empty()) {
    stack.resize (std::max (SIGSTKSZ, MINSIGSTKSZ) + 2*1024*1024);
    stack_t ss;
    ss.ss_sp = stack.data();
    ss.ss_flags = 0;
    ss.ss_size = stack.size();
    sigaltstack (&ss, nullptr);
  }
}


thread_local std::vector<uint8_t> CoreDumpSvc::s_stack;
