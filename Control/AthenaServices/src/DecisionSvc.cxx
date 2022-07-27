/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "DecisionSvc.h"

// Gaudi include files
#include "GaudiKernel/EventContext.h"
#include "PersistentDataModel/AthenaAttributeList.h"
#include <algorithm>

DecisionSvc::DecisionSvc(const std::string& name,
		       ISvcLocator* pSvcLocator ) :
  base_class(name, pSvcLocator),
  m_evtStore("StoreGateSvc",name),
#ifdef SIMULATIONBASE
  m_cutflowSvc("",name),
#else
  m_cutflowSvc("CutFlowSvc/CutFlowSvc",name),
#endif
  m_algstateSvc("AlgExecStateSvc",name)
{
}

DecisionSvc::~DecisionSvc()
{
}

StatusCode
DecisionSvc::initialize()
{
  // Decode the accept, required and veto Algorithms.
  // The logic is the following:
  //  a. The event is accepted if all lists are empty.
  //  b. The event is provisionally accepted if any Algorithm in the
  //     accept list
  //     has been executed and has indicated that its filter is passed. This
  //     provisional acceptance can be overridden by the other lists.
  //  c. The event is rejected unless all Algorithms in the required list have
  //     been executed and have indicated that their filter passed.
  //  d. The event is rejected if any Algorithm in the veto list has been
  //     executed and has indicated that its filter has passed.

  // Set to be listener for BeginRun
  // Declares top filters to ICutFlowSvc if (and only if) needed (David Cote, Sep 2010)

  // Get handle to exec state service for retrieving decisions
  ATH_CHECK( m_algstateSvc.retrieve() );

  ATH_CHECK(m_evtStore.retrieve());

  return StatusCode::SUCCESS;
}


StatusCode
DecisionSvc::finalize()
{
  ATH_MSG_INFO ("Finalized successfully.");
  return StatusCode::SUCCESS;
}

StatusCode
DecisionSvc::interpretAlgMap()
{
  StatusCode status = StatusCode::SUCCESS;

  // Fill up m_streamNames vector
  std::set<std::string> uniStreams;
  m_streamNames.clear();
  // first take list of streams from accept streams
  for(auto iter = m_stream_accept.begin();
           iter != m_stream_accept.end(); ++iter){
    uniStreams.insert(iter->first);
  }
  // then require
  for(auto iter = m_stream_require.begin();
           iter != m_stream_require.end(); ++iter){
    uniStreams.insert(iter->first);
  }
  // then veto
  for(auto iter = m_stream_veto.begin();
           iter != m_stream_veto.end(); ++iter){
    uniStreams.insert(iter->first);
  }
  std::copy(uniStreams.begin(),uniStreams.end(),std::back_inserter(m_streamNames));

  return status;
}


///////////////////////////////////////////////////////////////////
// Non-const methods:
///////////////////////////////////////////////////////////////////


StatusCode
DecisionSvc::addStream(const std::string& stream)
{
  StatusCode status = StatusCode::SUCCESS;
  if(m_frozen != true){
    // check if this stream already exist
    auto it = m_stream_accept.find(stream);
    if(it != m_stream_accept.end()){
      // ok, it exists, then do nothing
      ATH_MSG_WARNING("Stream name : " << stream << " already been registered!");
      status = StatusCode::FAILURE;
    }else{
      //if the stream doesn't exist yet, then insert it to the accept list with an empty vector of Algs
      std::vector<std::string> tmpvec;
      tmpvec.clear();
      ATH_MSG_INFO("Inserting stream: "<< stream  << " with no Algs");
      m_stream_accept.insert(std::make_pair(stream, tmpvec));
      status = StatusCode::SUCCESS;
    }
  }
  return status;
}

StatusCode
DecisionSvc::fillMap(std::map<std::string, std::vector<std::string> >& streamsModeMap,
                     const std::string& name, const std::string& stream)
{
  // check if this stream already exist
  std::map<std::string, std::vector<std::string> >::iterator it = streamsModeMap.find(stream);
  if(it != streamsModeMap.end()){
    // ok, it exists, then check if the algname was already been inserted


    bool algexist = false;
    // Check if alg already registered for this stream
    for(auto vit = (it->second).begin();
             vit != (it->second).end(); ++vit) {
      if((*vit) == name){
        algexist = true;
        // it seems the alg was already inserted, warn the user
        ATH_MSG_ERROR("Alg name : " << name
                       << " of stream " << stream
                       << " has already been registered!");
        return StatusCode::FAILURE;
      }
    }


    // So, if the stream exist but the alg has not been registered
    // update its content std::vector with a alg
    if(algexist == false){
      std::vector<std::string> &tmpvec = it->second;
      tmpvec.push_back(name);
    }

    //if the stream doesn't exist yet, then insert it
  } else {
    std::vector<std::string> tmpvec;
    tmpvec.push_back(name);
    streamsModeMap[stream] = std::move(tmpvec);
  }

  return StatusCode::SUCCESS;
}

StatusCode
DecisionSvc::addAcceptAlg(const std::string& name, const std::string& stream)
{
  // Add alg to list of streams
  auto it = std::find(m_streamNames.begin(), m_streamNames.end(),stream);
  if (it != m_streamNames.end()) m_streamNames.push_back(stream);
  // Fill map of stream to alg
  return fillMap(m_stream_accept,name,stream);
}

StatusCode
DecisionSvc::addRequireAlg(const std::string& name,
                           const std::string& stream)
{
  // Add alg to list of streams
  auto it = std::find(m_streamNames.begin(), m_streamNames.end(),stream);
  if (it != m_streamNames.end()) m_streamNames.push_back(stream);
  // Fill map of stream to alg
  return fillMap(m_stream_require,name,stream);
}

StatusCode
DecisionSvc::addVetoAlg(const std::string& name,
                        const std::string& stream)
{
  // Add alg to list of streams
  auto it = std::find(m_streamNames.begin(), m_streamNames.end(),stream);
  if (it != m_streamNames.end()) m_streamNames.push_back(stream);
  // Fill map of stream to alg
  return fillMap(m_stream_veto,name,stream);
}

// Utility interface to retrieve unique list of streams
const std::vector<std::string> DecisionSvc::getStreams() const
{
  return m_streamNames;
}

// Utility interface 
const std::vector<std::string>
DecisionSvc::getAcceptAlgs(const std::string& stream) const
{
  auto iter = m_stream_accept.find(stream);
  if (iter != m_stream_accept.end()) return iter->second;
  return std::vector<std::string>();
}

// Utility interface
const std::vector<std::string>
DecisionSvc::getRequireAlgs(const std::string& stream) const
{
  auto iter = m_stream_require.find(stream);
  if (iter != m_stream_require.end()) return iter->second;
  return std::vector<std::string>();
}

// Utility interface
const std::vector<std::string>
DecisionSvc::getVetoAlgs(const std::string& stream) const
{
  auto iter = m_stream_veto.find(stream);
  if (iter != m_stream_veto.end()) return iter->second;
  return std::vector<std::string>();
}



///////////////////////////////////////////////////////////////////
// Const methods:
///////////////////////////////////////////////////////////////////


bool
DecisionSvc::isEventAccepted( const std::string& stream ) const
{
  EventContext ec;
  ec.setValid(false);
  ec.setSlot(0);
  return isEventAccepted(stream,ec);
}

bool
DecisionSvc::isEventAccepted( const std::string& stream, 
                              const EventContext& ectx) const
{

  ATH_MSG_DEBUG("In DecisionSvc::isEventAccepted( " << stream << " )");

  // By construction a stream is accepted
  bool result = true;
  
  bool found_accept = false;
  bool found_require = false;
  bool found_veto = false;

  //Loop over all streams of accept type and find the one of interest
  auto itAlgs = m_stream_accept.find(stream);
  if(itAlgs != m_stream_accept.end()){
    found_accept = true;
    // get a handle of the streams' algos vector
    const auto &vecAlgs = itAlgs->second;
    // Loop over all Algorithms in the list to see
    // whether any have been executed and have their filter
    // passed flag set. Any match causes the event to be
    // provisionally accepted.
    if ( ! vecAlgs.empty( ) ) {
      result = false;
      for (auto it = vecAlgs.begin(); it != vecAlgs.end(); it++) {
        bool isE,fp;
        isE = (m_algstateSvc->algExecState(*it,ectx).state() == AlgExecState::State::Done);
        fp = m_algstateSvc->algExecState(*it,ectx).filterPassed();
        if (isE && fp) {
	  result = true;
	  break;
	}
      }
    }
  }

  //Loop over all streams of require type and find the one of interest
  itAlgs = m_stream_require.find(stream);
  if(itAlgs != m_stream_require.end()){
    found_require = true;
    // get a handle of the streams' algos vector
    const auto &vecAlgs = itAlgs->second;
    // Loop over all Algorithms in the list to see
    // whether any have been executed and have their filter
    // passed flag set. Any match causes the event to be
    // provisionally accepted.
    if ( ! vecAlgs.empty( ) ) {
      for (auto it = vecAlgs.begin(); it != vecAlgs.end(); it++) {
        bool isE,fp;
        isE = (m_algstateSvc->algExecState(*it,ectx).state() == AlgExecState::State::Done);
        fp = m_algstateSvc->algExecState(*it,ectx).filterPassed();
        if (!isE || !fp) {
          result = false;
          break;
        }
      }
    }
  }

  //Loop over all streams of veto type and find the one of interest
  itAlgs = m_stream_veto.find(stream);
  if(itAlgs != m_stream_veto.end()){
    found_veto = true;
    // get a handle of the streams' algos vector
    const auto &vecAlgs = itAlgs->second;
    // Loop over all Algorithms in the list to see
    // whether any have been executed and have their filter
    // passed flag set. Any match causes the event to be
    // provisionally accepted.
    if ( ! vecAlgs.empty( ) ) {
      for (auto it = vecAlgs.begin(); it != vecAlgs.end(); it++) {
        bool isE,fp;
        isE = (m_algstateSvc->algExecState(*it,ectx).state() == AlgExecState::State::Done);
        fp = m_algstateSvc->algExecState(*it,ectx).filterPassed();
	if ( isE && fp ) {
          result = false;
          break;
	}
      }
    }
  }

  if(found_accept == false && found_require == false && found_veto == false){
    ATH_MSG_DEBUG("Stream: " << stream << " not found registered in DecisionSvc -- accepting event by default ");
  }

  return result;
}

StatusCode DecisionSvc::start()
{
  ATH_MSG_DEBUG("in start");
  CHECK( this->interpretAlgMap() );
  m_frozen = true;

  // lambda to return true if second element is non-empty
  auto teststreams = [](const auto& m)
    {
      for (const auto& p : m)
	if (!p.second.empty()) return true;
      return false;
    };
  
  //Retrieve CutFlowSvc if (and only if) needed
  if( teststreams(m_stream_accept) || teststreams(m_stream_require) || teststreams(m_stream_veto)) {
    if (!m_cutflowSvc.empty())
      {
        if (m_cutflowSvc.retrieve().isFailure())
          {
            ATH_MSG_ERROR("Cannot get ICutFlowSvc interface.");
            return StatusCode::RECOVERABLE;
          }
      }
  }

  //Now that everything is said and done, match filters with stream and logic in CutFlowSvc
  StatusCode sc ATLAS_THREAD_SAFE = DeclareToCutFlowSvc();
  //            ^ FIXME: Should really mark start() as ATLAS_NOT_THREAD_SAFE but that is currently
  //              not possible with the thread-checker.

  return sc;
}

StatusCode DecisionSvc::DeclareToCutFlowSvc ATLAS_NOT_THREAD_SAFE ()
{
  // Declare all decisions to CutFlowSvc for bookkeeping
  
  // Loop over all streams
  for(auto ait  = m_streamNames.begin();
           ait != m_streamNames.end(); ++ait) {
    const std::string &streamName=*ait;
    const std::vector<std::string> *accFilt=nullptr;
    const std::vector<std::string> *reqFilt=nullptr;
    const std::vector<std::string> *vetFilt=nullptr;
    std::vector< const std::vector<std::string>* > totFilt;
    if (auto itr = m_stream_accept.find(streamName); itr != m_stream_accept.end())  {
      accFilt = &itr->second;
      totFilt.push_back(accFilt);
    }
    if (auto itr = m_stream_require.find(streamName); itr != m_stream_require.end()) {
      reqFilt = &itr->second;
      totFilt.push_back(reqFilt);
    }
    if (auto itr = m_stream_veto.find(streamName); itr != m_stream_veto.end())    {
      vetFilt = &itr->second;
      totFilt.push_back(vetFilt);
    }

    // Now build logicalKey as string of filt ||, &&, ! based on 
    // whether it is accept, require, veto
    std::string accstring(""), reqstring(""), vetstring("");
    if(accFilt){
       for (auto it = accFilt->begin(); it != accFilt->end(); ++it) {
         if(!accstring.empty()) accstring += "||";
         else accstring += '(';
         accstring+=*it;
       }
    }
    if(reqFilt){
       for (auto it = reqFilt->begin(); it != reqFilt->end(); ++it) {
         if(!reqstring.empty()) reqstring += "&&";
         else reqstring += '(';
         reqstring+=*it;
       }
    }
    if(vetFilt){
       for (auto it = vetFilt->begin(); it != vetFilt->end(); ++it) {
         if(!vetstring.empty()) vetstring += "||";
         else vetstring += '(';
         vetstring+=*it;
       }
    }
    std::string logicalKey("");
    if(!accstring.empty()) {
      accstring += ')';
      logicalKey += accstring;
    }
    if(!reqstring.empty()) {
      reqstring += ')';
      if (!logicalKey.empty()) logicalKey += "&&";
      logicalKey += reqstring;
    }
    if(!vetstring.empty()) {
      vetstring += ')';
      if (!logicalKey.empty()) logicalKey += "&&";
      logicalKey += '!';
      logicalKey += vetstring;
    }
    // If no filters, mark as PasThru
    if (logicalKey.empty()) logicalKey="PassThru";
    ATH_MSG_DEBUG("stream " << streamName << " uses logic " << logicalKey);

    // Now actually declare to the cutflowsvc
    for(auto vec  = totFilt.begin();
             vec != totFilt.end(); ++vec) {
      ATH_MSG_DEBUG("Declaring logic " << logicalKey << " for " << streamName);
      for (auto filter  = (*vec)->begin();
                filter != (*vec)->end(); ++filter) {
        if(!m_cutflowSvc.empty()) {m_cutflowSvc->registerTopFilter( (*filter), logicalKey, 2, streamName, true ); } // TODO: validate
      }
    }
  }

  return StatusCode::SUCCESS;
}
