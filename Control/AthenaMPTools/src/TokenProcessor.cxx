/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "TokenProcessor.h"
#include "copy_file_icc_hack.h"
#include "AthenaInterprocess/ProcessGroup.h"

#include "AthenaKernel/IEventShare.h"
#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/IIoComponentMgr.h"
#include "GaudiKernel/IFileMgr.h"
#include "GaudiKernel/IChronoStatSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/FileIncident.h"
#include "GaudiKernel/Timing.h"

#include <sys/stat.h>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <stdint.h>
#include <stdexcept>
#include <queue>

#include "yampl/SocketFactory.h"

struct ShareEventHeader {
  long evtSeqNumber;
  long fileSeqNumber;
  std::size_t evtSize;
  std::size_t evtOffset;
  unsigned int evtCoreStatusFlag;
  uint32_t evtTerm1;
  uint32_t evtTerm2;
};


TokenProcessor::TokenProcessor(const std::string& type
			       , const std::string& name
			       , const IInterface* parent)
  : AthenaMPToolBase(type,name,parent)
  , m_isPileup(false)
  , m_rankId(-1)
  , m_nEventsBeforeFork(0)
  , m_chronoStatSvc("ChronoStatSvc", name)
  , m_evtShare(0)
  , m_channel2Scatterer("")
  , m_channel2EvtSel("")
  , m_sharedRankQueue(0)
  , m_socketFactory(0)
  , m_socket2Scatterer(0)
{
  declareInterface<IAthenaMPTool>(this);

  declareProperty("IsPileup",m_isPileup);
  declareProperty("EventsBeforeFork",m_nEventsBeforeFork);
  declareProperty("Channel2Scatterer", m_channel2Scatterer);
  declareProperty("Channel2EvtSel", m_channel2EvtSel);

  m_subprocDirPrefix = "worker_";
}

TokenProcessor::~TokenProcessor()
{
}

StatusCode TokenProcessor::initialize()
{
  msg(MSG::DEBUG) << "In initialize" << endreq;
  if(m_isPileup) {
    m_evtProcessor = ServiceHandle<IEventProcessor>("PileUpEventLoopMgr",name());
    msg(MSG::INFO) << "The job running in pileup mode" << endreq;
  }
  else {
    msg(MSG::INFO) << "The job running in non-pileup mode" << endreq;
  }

  StatusCode sc = AthenaMPToolBase::initialize();
  if(!sc.isSuccess())
    return sc;

  sc = serviceLocator()->service(m_evtSelName,m_evtShare);
  if(sc.isFailure() || m_evtShare==0) {
    msg(MSG::ERROR) << "Error retrieving IEventShare" << endreq;
    return StatusCode::FAILURE;
  }
  
  sc = m_chronoStatSvc.retrieve();
  if (!sc.isSuccess()) {
    msg(MSG::ERROR) << "Cannot get ChronoStatSvc." << endreq;
    return StatusCode::FAILURE;
  }
  
  return StatusCode::SUCCESS;
}

StatusCode TokenProcessor::finalize()
{
  delete m_sharedRankQueue;
  return StatusCode::SUCCESS;
}

int TokenProcessor::makePool(int, int nprocs, const std::string& topdir)
{
  msg(MSG::DEBUG) << "In makePool " << getpid() << endreq;

  if(nprocs==0 || nprocs<-1) {
    msg(MSG::ERROR) << "Invalid value for the nprocs parameter: " << nprocs << endreq;
    return -1;
  }

  if(topdir.empty()) {
    msg(MSG::ERROR) << "Empty name for the top directory!" << endreq;
    return -1;
  }

  m_nprocs = (nprocs==-1?sysconf(_SC_NPROCESSORS_ONLN):nprocs);
  m_subprocTopDir = topdir;

  // Create rank queue and fill it
  std::ostringstream rankQueueName;
  rankQueueName << "TokenProcessor_RankQueue_" << getpid();
  m_sharedRankQueue = new AthenaInterprocess::SharedQueue(rankQueueName.str(),m_nprocs,sizeof(int));
  for(int i=0; i<m_nprocs; ++i)
    if(!m_sharedRankQueue->send_basic<int>(i)) {
      msg(MSG::ERROR) << "Unable to send int to the ranks queue!" << endreq;
      return -1;
    }

  // Create the process group and map_async bootstrap
  m_processGroup = new AthenaInterprocess::ProcessGroup(m_nprocs);
  msg(MSG::INFO) << "Created Pool of " << m_nprocs << " worker processes" << endreq;
  if(mapAsyncFlag(AthenaMPToolBase::FUNC_BOOTSTRAP))
    return -1;
  msg(MSG::INFO) << "Workers bootstraped" << endreq; 

  return m_nprocs;
}

StatusCode TokenProcessor::exec()
{
  msg(MSG::DEBUG) << "In exec " << getpid() << endreq;

  if(mapAsyncFlag(AthenaMPToolBase::FUNC_EXEC))
    return StatusCode::FAILURE;
  msg(MSG::INFO) << "Workers started processing events" << endreq;

  return StatusCode::SUCCESS;
}

StatusCode TokenProcessor::wait_once(int& numFinishedProc)
{
  StatusCode sc = AthenaMPToolBase::wait_once(numFinishedProc);
  AthenaInterprocess::ProcessResult* presult(0);
  if(sc.isFailure()) {
    // Try to start a new process
    if(startProcess().isSuccess()) {
      msg(MSG::INFO) << "Successfully started new process" << endreq;
      numFinishedProc--;
    }
    else
      msg(MSG::WARNING) << "Failed to start new process" << endreq;
    
    if(m_processGroup->getChildren().size()) {
      msg(MSG::INFO) << "The process group continues with " << m_processGroup->getChildren().size() << " processes" << endreq;
      return StatusCode::SUCCESS;
    }
    msg(MSG::ERROR) << "No more processes in the group!" << endreq;
  }
  else {
    // Pull one result and decode it if necessary
    presult = m_processGroup->pullOneResult();
    int res(0);
    if(presult && (unsigned)(presult->output.size)>sizeof(int))
      res = decodeProcessResult(presult,true);
    if(presult) free(presult->output.data);
    delete presult;
    if(res) return StatusCode::FAILURE;
  }
  return sc;
}

void TokenProcessor::reportSubprocessStatuses()
{
  msg(MSG::INFO) << "Statuses of event processors" << endreq;
  const std::vector<AthenaInterprocess::ProcessStatus>& statuses = m_processGroup->getStatuses();
  for(size_t i=0; i<statuses.size(); ++i) {
    // Get the number of events processed by this worker
    std::map<pid_t,int>::const_iterator it = m_nProcessedEvents.find(statuses[i].pid);
    msg(MSG::INFO) << "*** Process PID=" << statuses[i].pid 
		   << ". Status " << ((statuses[i].exitcode)?"FAILURE":"SUCCESS") 
		   << ". Number of events processed: ";
    if(it==m_nProcessedEvents.end())
      msg(MSG::INFO) << "N/A" << endreq;
    else
      msg(MSG::INFO) << it->second << endreq;
  }
}

void TokenProcessor::subProcessLogs(std::vector<std::string>& filenames)
{
  filenames.clear();
  for(int i=0; i<m_nprocs; ++i) {
    std::ostringstream workerIndex;
    workerIndex << i;
    boost::filesystem::path worker_rundir(m_subprocTopDir);
    worker_rundir /= boost::filesystem::path(m_subprocDirPrefix+workerIndex.str());
    filenames.push_back(worker_rundir.string()+std::string("/AthenaMP.log"));
  }
}

AthenaMP::AllWorkerOutputs_ptr TokenProcessor::generateOutputReport()
{
  AthenaMP::AllWorkerOutputs_ptr jobOutputs(new AthenaMP::AllWorkerOutputs());
  return jobOutputs;
}

AthenaInterprocess::ScheduledWork* TokenProcessor::bootstrap_func()
{
  int* errcode = new int(1); // For now use 0 success, 1 failure
  AthenaInterprocess::ScheduledWork* outwork = new AthenaInterprocess::ScheduledWork;
  outwork->data = (void*)errcode;
  outwork->size = sizeof(int);
  // ...
  // (possible) TODO: extend outwork with some error message, which will be eventually
  // reported in the master proces
  // ...

  // ________________________ Get RankID ________________________
  //
  if(!m_sharedRankQueue->receive_basic<int>(m_rankId)) {
    msg(MSG::ERROR) << "Unable to get rank ID!" << endreq;
    return outwork;
  }
  std::ostringstream workindex;
  workindex<<m_rankId;

  // ________________________ Worker dir: mkdir ________________________
  boost::filesystem::path worker_rundir(m_subprocTopDir);
  worker_rundir /= boost::filesystem::path(m_subprocDirPrefix+workindex.str());
  // TODO: this "worker_" can be made configurable too

  if(mkdir(worker_rundir.string().c_str(),S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH)==-1) {
    msg(MSG::ERROR) << "Unable to make worker run directory: " << worker_rundir.string() << ". " << strerror(errno) << endreq;
    return outwork;
  }

  // ________________________ Redirect logs ________________________
  if(redirectLog(worker_rundir.string()))
    return outwork;

  msg(MSG::INFO) << "Logs redirected in the AthenaMP event worker PID=" << getpid() << endreq;

  // ________________________ Update Io Registry ____________________________
  if(updateIoReg(worker_rundir.string()))
    return outwork;

  msg(MSG::INFO) << "Io registry updated in the AthenaMP event worker PID=" << getpid() << endreq;

  // ________________________ SimParams & DigiParams ____________________________
  boost::filesystem::path abs_worker_rundir = boost::filesystem::absolute(worker_rundir);
  if(boost::filesystem::is_regular_file("SimParams.db"))
    COPY_FILE_HACK("SimParams.db", abs_worker_rundir.string()+"/SimParams.db");
  if(boost::filesystem::is_regular_file("DigitParams.db"))
    COPY_FILE_HACK("DigitParams.db", abs_worker_rundir.string()+"/DigitParams.db");

  // _______________________ Handle saved PFC (if any) ______________________
  if(handleSavedPfc(abs_worker_rundir))
    return outwork;

  // ________________________  reopen descriptors ____________________________
  if(reopenFds())
    return outwork;

  msg(MSG::INFO) << "File descriptors re-opened in the AthenaMP event worker PID=" << getpid() << endreq;

  
  // ________________________ Make Shared RAW Reader Client ________________________
  if(!m_evtShare->makeClient(m_rankId).isSuccess()) {
    msg(MSG::ERROR) << "Failed to make the event selector a share client" << endreq;
    return outwork;
  }
  else {
    msg(MSG::DEBUG) << "Successfully made the event selector a share client" << endreq;
  }

  // ________________________ I/O reinit ________________________
  if(!m_ioMgr->io_reinitialize().isSuccess()) {
    msg(MSG::ERROR) << "Failed to reinitialize I/O" << endreq;
    return outwork;
  } else {
    msg(MSG::DEBUG) << "Successfully reinitialized I/O" << endreq;
  }

  // ________________________ Event selector restart ________________________
  IService* evtSelSvc = dynamic_cast<IService*>(m_evtSelector);
  if(!evtSelSvc) {
    msg(MSG::ERROR) << "Failed to dyncast event selector to IService" << endreq;
    return outwork;
  }
  if(!evtSelSvc->start().isSuccess()) {
    msg(MSG::ERROR) << "Failed to restart the event selector" << endreq;
    return outwork;
  } else {
    msg(MSG::DEBUG) << "Successfully restarted the event selector" << endreq;
  }

  // ________________________ Restart background event selectors in pileup jobs ________________________
  if(m_isPileup) {
    const std::list<IService*>& service_list = serviceLocator()->getServices();
    std::list<IService*>::const_iterator itSvc = service_list.begin(),
      itSvcLast = service_list.end();
    for(;itSvc!=itSvcLast;++itSvc) {
      IEvtSelector* evtsel = dynamic_cast<IEvtSelector*>(*itSvc);
      if(evtsel && (evtsel != m_evtSelector)) {
	if((*itSvc)->start().isSuccess())
	  msg(MSG::DEBUG) << "Restarted event selector " << (*itSvc)->name() << endreq;
	else {
	  msg(MSG::ERROR) << "Failed to restart event selector " << (*itSvc)->name() << endreq;
	  return outwork;
	}
      }
    }
  }

  // ________________________ Worker dir: chdir ________________________
  if(chdir(worker_rundir.string().c_str())==-1) {
    msg(MSG::ERROR) << "Failed to chdir to " << worker_rundir.string() << endreq;
    return outwork;
  }

  // Declare success and return
  *errcode = 0;
  return outwork;
}

AthenaInterprocess::ScheduledWork* TokenProcessor::exec_func()
{
  msg(MSG::INFO) << "Exec function in the AthenaMP worker PID=" << getpid() << endreq;

  int nEvt(1);
  int nEventsProcessed(0);

  bool all_ok(true);
  long eventNumber(0);
  std::queue<std::string> queueTokens;

  // Get the yampl connection channels
  m_socketFactory = new yampl::SocketFactory();
  m_socket2Scatterer = m_socketFactory->createClientSocket(yampl::Channel(m_channel2Scatterer.value(),yampl::LOCAL_PIPE),yampl::MOVE_DATA);
  msg(MSG::DEBUG) << "Created CLIENT socket to the Scatterer: " << m_channel2Scatterer.value() << endreq;
  std::ostringstream pidstr;
  pidstr << getpid();
  std::string socket2EvtSelName = m_channel2EvtSel.value() + std::string("_") + pidstr.str();
  yampl::ISocket* socket2EvtSel = m_socketFactory->createClientSocket(yampl::Channel(socket2EvtSelName,yampl::LOCAL_PIPE),yampl::COPY_DATA);
  msg(MSG::DEBUG) << "Created CLIENT socket to the Tool: " << socket2EvtSelName << endreq;

  // Get the IncidentSvc
  IIncidentSvc* p_incidentSvc(0);
  if(!serviceLocator()->service("IncidentSvc", p_incidentSvc).isSuccess()) {
    msg(MSG::ERROR) << "Unable to retrieve IncidentSvc" << endreq;
    all_ok = false;
  }

  // Construct a "welcome" message to be sent to the TokenScatterer
  std::string ping = pidstr.str() + std::string(" ready for event processing");
  void* message2scatterer = malloc(ping.size());
  memcpy(message2scatterer,ping.data(),ping.size());
  m_socket2Scatterer->send(message2scatterer,ping.size());
  msg(MSG::DEBUG) << "Sent a welcome message to the Scatterer" << endreq;

  while(all_ok) {
    // Get the response - list of tokens - from the scatterer. 
    // The format of the response: | ResponseSize | RangeID, | evtToken[,evtToken] |
    char *responseBuffer(0);
    ssize_t responseSize = m_socket2Scatterer->recv(responseBuffer);
    // If response size is 0 then break the loop
    if(responseSize==1) {
      msg(MSG::DEBUG) << "Empty range received. Terminating the loop" << endreq;
      break;
    }

    std::string responseStr(responseBuffer,responseSize);
    msg(MSG::DEBUG) << "Received response from the Scatterer : " << responseStr << endreq;

    // Start timing
    System::ProcessTime time_start = System::getProcessTime();

    size_t startpos(0);
    size_t endpos = responseStr.find(',');
    while(endpos!=std::string::npos) {
      queueTokens.push(responseStr.substr(startpos,endpos-startpos));
      startpos = endpos+1;
      endpos = responseStr.find(',',startpos);
    }
    queueTokens.push(responseStr.substr(startpos));
    // Actually the first element in the tokens queue is the RangeID. Get it
    std::string rangeID = queueTokens.front();
    queueTokens.pop();
    msg(MSG::DEBUG) << "Received RangeID=" << rangeID << endreq;
    // Fire an incident
    if(!queueTokens.empty()) {
      p_incidentSvc->fireIncident(FileIncident(name(),"NextEventRange",rangeID));

      // Time to report the previous output
      if(!m_outputFileReport.empty()) {
	message2scatterer = malloc(m_outputFileReport.size());
	memcpy(message2scatterer,m_outputFileReport.data(),m_outputFileReport.size());
	m_socket2Scatterer->send(message2scatterer,m_outputFileReport.size());
	m_outputFileReport.clear();
      }
    }

    // Pass the tokens over to the Event Selector one at a time
    while(!queueTokens.empty()) {
      std::string strToken = queueTokens.front();
      msg(MSG::DEBUG) << "Processing the Token : " << strToken << endreq;
      queueTokens.pop();

      std::size_t num = strToken.size();

      ShareEventHeader evtH;
      evtH.evtSeqNumber = eventNumber++;
      evtH.fileSeqNumber = 0;
      evtH.evtSize = num;
      evtH.evtOffset = 0;
      evtH.evtCoreStatusFlag = 0; // ???
      evtH.evtTerm1 = *(static_cast<const uint32_t*>((const void*)strToken.data()) + num / sizeof(uint32_t) - 1);
      evtH.evtTerm2 = *(static_cast<const uint32_t*>((const void*)strToken.data()) + num / sizeof(uint32_t) - 2);

      int messageSize = num+sizeof(evtH);
      void* message = malloc(messageSize);
      memcpy(message,(void*)&evtH,sizeof(evtH));
      memcpy((char*)message+sizeof(evtH),strToken.data(),num);

      socket2EvtSel->send(message,num+sizeof(evtH));
      msg(MSG::DEBUG) << "Message sent to the event selector" << endreq;

      // Process the event
      m_chronoStatSvc->chronoStart("AthenaMP_nextEvent");
      StatusCode sc = m_evtProcessor->nextEvent(nEvt++);
      nEventsProcessed++;
      if(sc.isFailure()){
	msg(MSG::ERROR) << "Unable to process the event" << endreq;
	all_ok=false;
      }
      else
	msg(MSG::DEBUG)<< "Event processed" << endreq;
      m_chronoStatSvc->chronoStop("AthenaMP_nextEvent"); 

      // Get the response posted by the event selector
      char *pong(0); // can be something else
      socket2EvtSel->recv(pong);
      msg(MSG::DEBUG) << "Event selector response received from the channel" << endreq;

      if(!all_ok) break;
    } // Loop over tokens in the event range

    std::string strOutpFile;
    if(all_ok) {
      // Get the full path of the event range output file
      for(boost::filesystem::directory_iterator fdIt(boost::filesystem::current_path()); fdIt!=boost::filesystem::directory_iterator(); fdIt++) {
	if(fdIt->path().string().find(rangeID)!=std::string::npos) {
	  if(!strOutpFile.empty()) {
	    msg(MSG::ERROR) << "More than one file containing RangeID=" << rangeID << " found in the run dir" << endreq;
	    all_ok = false;
	    break;
	  }
	  strOutpFile = fdIt->path().string();
	}
      }
    }

    if(!all_ok) break;

    // Stop timing
    System::ProcessTime time_delta = System::getProcessTime() - time_start;
    
    // Prepare the output report
    if(!strOutpFile.empty()) {
      // We need to combine the output file name with
      // 1. RangeID (requested by JEDI)
      // 2. CPU time
      // 3. Wall time
      std::ostringstream outputReportStream;
      outputReportStream << strOutpFile << "," << rangeID 
			 << ",CPU:" << time_delta.cpuTime<System::Sec>()
			 << ",WALL:" << time_delta.elapsedTime<System::Sec>();
      m_outputFileReport = outputReportStream.str();
    }

    // Request the next available range
    message2scatterer = malloc(ping.size());
    memcpy(message2scatterer,ping.data(),ping.size());
    m_socket2Scatterer->send(message2scatterer,ping.size());
    msg(MSG::DEBUG) << "Sent a message to the scatterer: " << ping << endreq;
  } // Main "event loop"

  if(all_ok) {
    if(m_evtProcessor->executeRun(0).isFailure()) {
      msg(MSG::ERROR) << "Could not finalize the Run" << endreq;
      all_ok=false;
    }
  }

  int errcode = (all_ok?0:1); // For now use 0 success, 1 failure
  AthenaMPToolBase::Func_Flag func = AthenaMPToolBase::FUNC_EXEC;
  // Return value: "ERRCODE|Func_Flag|NEvt"
  int outsize = 2*sizeof(int)+sizeof(AthenaMPToolBase::Func_Flag);
  void* outdata = malloc(outsize);
  memcpy(outdata,&errcode,sizeof(int));
  memcpy((char*)outdata+sizeof(int),&func,sizeof(func));
  memcpy((char*)outdata+sizeof(int)+sizeof(func),&nEventsProcessed,sizeof(int));
  AthenaInterprocess::ScheduledWork* outwork = new AthenaInterprocess::ScheduledWork;
  outwork->data = outdata;
  outwork->size = outsize;
  // ...
  // (possible) TODO: extend outwork with some error message, which will be eventually
  // reported in the master proces
  // ...

  delete socket2EvtSel;

  return outwork;
}

AthenaInterprocess::ScheduledWork* TokenProcessor::fin_func()
{
  msg(MSG::INFO) << "Fin function in the AthenaMP worker PID=" << getpid() << endreq;

  // We are not able to use private data members after the appMgr has been finalized
  yampl::ISocket* socket2Scatterer(m_socket2Scatterer);
  yampl::ISocketFactory* socketFactory(m_socketFactory);
  std::string outputFileReport(m_outputFileReport);

  bool all_ok(true);

  if(m_appMgr->stop().isFailure()) {
    msg(MSG::ERROR) << "Unable to stop AppMgr" << endreq; 
    all_ok=false;
  }
  else { 
    if(m_appMgr->finalize().isFailure()) {
      msg(MSG::ERROR) << "Unable to finalize AppMgr" << endreq;
      all_ok=false;
    }
  }

  // Report the last output file
  if(!outputFileReport.empty()) {
    void* message2scatterer = malloc(outputFileReport.size());
    memcpy(message2scatterer,outputFileReport.data(),outputFileReport.size());
    socket2Scatterer->send(message2scatterer,outputFileReport.size());
  }

  int errcode = (all_ok?0:1); // For now use 0 success, 1 failure
  int nEvt = -1;
  AthenaMPToolBase::Func_Flag func = AthenaMPToolBase::FUNC_FIN;
  // Return value: "ERRCODE|Func_Flag|NEvt"  (Here NEvt=-1)
  int outsize = 2*sizeof(int)+sizeof(AthenaMPToolBase::Func_Flag);
  void* outdata = malloc(outsize);
  memcpy(outdata,&errcode,sizeof(int));
  memcpy((char*)outdata+sizeof(int),&func,sizeof(func));
  memcpy((char*)outdata+sizeof(int)+sizeof(func),&nEvt,sizeof(int));
  AthenaInterprocess::ScheduledWork* outwork = new AthenaInterprocess::ScheduledWork;
  outwork->data = outdata;
  outwork->size = outsize;

  delete socket2Scatterer;
  delete socketFactory;

  return outwork;
}

int TokenProcessor::decodeProcessResult(const AthenaInterprocess::ProcessResult* presult, bool doFinalize)
{
  if(!presult) return 0;
  const AthenaInterprocess::ScheduledWork& output = presult->output;
  msg(MSG::DEBUG) << "Decoding the output of PID=" << presult->pid << " with the size=" << output.size << endreq;
  if(output.size!=2*sizeof(int)+sizeof(AthenaMPToolBase::Func_Flag)) return 0;
  
  AthenaMPToolBase::Func_Flag func;
  memcpy(&func,(char*)output.data+sizeof(int),sizeof(func));

  if(func==AthenaMPToolBase::FUNC_EXEC) {
    // Store the number of processed events
    int nevt(0);
    memcpy(&nevt,(char*)output.data+sizeof(int)+sizeof(func),sizeof(int));
    m_nProcessedEvents[presult->pid]=nevt;
    msg(MSG::DEBUG) << "PID=" << presult->pid << " processed " << nevt << " events" << endreq;

    if(doFinalize) {
      // Add PID to the finalization queue
      m_finQueue.push(presult->pid);
      msg(MSG::DEBUG) << "Added PID=" << presult->pid << " to the finalization queue" << endreq;

      // If this is the only element in the queue then start its finalization
      // Otherwise it has to wait its turn until all previous processes have been finalized
      if(m_finQueue.size()==1) {
        if(mapAsyncFlag(AthenaMPToolBase::FUNC_FIN,presult->pid)
           || m_processGroup->map_async(0,0,presult->pid)) {
          msg(MSG::ERROR) << "Problem scheduling finalization on PID=" << presult->pid << endreq;
          return 1;
        }
        else {
          msg(MSG::DEBUG) << "Scheduled finalization of PID=" << presult->pid << endreq;
        }
      }
    }
  }
  else if(doFinalize && func==AthenaMPToolBase::FUNC_FIN) {
    msg(MSG::DEBUG) << "Finished finalization of PID=" << presult->pid << endreq;
    pid_t pid = m_finQueue.front();
    if(pid==presult->pid) {
      // pid received as expected. Remove it from the queue
      m_finQueue.pop();
      msg(MSG::DEBUG) << "PID=" << presult->pid << " removed from the queue" << endreq;
      // Schedule finalization of the next processe in the queue
      if(m_finQueue.size()) {
        if(mapAsyncFlag(AthenaMPToolBase::FUNC_FIN,m_finQueue.front())
           || m_processGroup->map_async(0,0,m_finQueue.front())) {
          msg(MSG::ERROR) << "Problem scheduling finalization on PID=" << m_finQueue.front() << endreq;
          return 1;
        }
        else  {
          msg(MSG::DEBUG) << "Scheduled finalization of PID=" << m_finQueue.front() << endreq;
        }
      }
    }
    else {
      // Error: unexpected pid received from presult
      msg(MSG::ERROR) << "Finalized PID=" << presult->pid << " while PID=" << pid << " was expected" << endreq;
      return 1;
    }
  }
  return 0;
}

StatusCode TokenProcessor::startProcess()
{
  m_nprocs++;

  // Create a rank for the new process
  if(!m_sharedRankQueue->send_basic<int>(m_nprocs-1)) {
    msg(MSG::WARNING) << "Unable to send int to the ranks queue!" << endreq;
    return StatusCode::FAILURE;
  }
  
  pid_t pid = m_processGroup->launchProcess();
  if(pid==0) {
    msg(MSG::WARNING) << "Unable to start new process" << endreq;
    return StatusCode::FAILURE;
  }
  
  if(mapAsyncFlag(AthenaMPToolBase::FUNC_BOOTSTRAP,pid)
     || mapAsyncFlag(AthenaMPToolBase::FUNC_EXEC,pid)) {
    msg(MSG::WARNING) << "Unable to map work on the new process" << endreq;
    return StatusCode::FAILURE;
  }
  
  if(m_processGroup->map_async(0,0,pid)){
    msg(MSG::WARNING) << "Unable to set exit to the new process" << endreq;
    return StatusCode::FAILURE;
  }

 return StatusCode::SUCCESS;
}
