/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef IOVDBTESTALG_PIXELDEDXTESTALG_H
#define IOVDBTESTALG_PIXELDEDXTESTALG_H


#include "GaudiKernel/Algorithm.h"
#define NCLASS 5
#define NPAR   9

#include "StoreGate/DataHandle.h"
#include "GaudiKernel/IIncidentListener.h"
#include "AthenaKernel/IOVSvcDefs.h"

class StoreGateSvc;
class EventInfo;
class IIOVRegistrationSvc;
class IAthenaOutputStreamTool;

/**
 ** Algorithm to test writing conditions data and reading them back.
 **/

class PixeldEdxTestAlg: public Algorithm 
{
public:
    PixeldEdxTestAlg (const std::string& name, ISvcLocator* pSvcLocator);
    ~PixeldEdxTestAlg();

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:

    StatusCode createCondObjects();
    StatusCode printCondObjects();
    StatusCode streamOutCondObjects();
    StatusCode registerCondObjects();
    StatusCode readWithBeginRun();
	StatusCode testCallBack(  IOVSVC_CALLBACK_ARGS  );

    StringProperty            m_streamName;
    std::string               m_tagID;    
    StoreGateSvc*             m_sgSvc;
    StoreGateSvc*             m_detStore;
IIOVRegistrationSvc*      m_regSvc;
    IAthenaOutputStreamTool*  m_streamer;
  std::string m_filename;
};


#endif 
