///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// RootOutputStreamTool.cxx 
// Implementation file for class Athena::RootOutputStreamTool
// Author Peter van Gemmeren <gemmeren@anl.gov>
// Author: S.Binet<binet@cern.ch>
/////////////////////////////////////////////////////////////////// 

// AthenaRootComps includes
#include "RootOutputStreamTool.h"
#include "RootSvc.h"
#include "RootConnection.h"
#include "RootBranchAddress.h"

// stl
#include "CxxUtils/unordered_set.h" // move to stl

// Gaudi
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/IRegistry.h"
#include "GaudiKernel/DataObject.h"

// Athena
#include "AthenaKernel/IClassIDSvc.h"
#include "AthenaRootKernel/IIoSvc.h"
#include "StoreGate/StoreGateSvc.h"
#include "SGTools/BuiltinsClids.h"

namespace Athena {

RootOutputStreamTool::RootOutputStreamTool(const std::string& type,
                                           const std::string& name,
                                           const IInterface* parent) : 
  ::AthAlgTool(type, name, parent),
  m_storeSvc("StoreGateSvc", name),
  m_conversionSvc("Athena::RootCnvSvc/AthenaRootCnvSvc", name),
  m_clidSvc("ClassIDSvc", name) 
{
  // Declare IAthenaOutputStreamTool interface
  declareInterface<IAthenaOutputStreamTool>(this);
  // Properties
  declareProperty("Store", 
                  m_storeSvc,
                  "Store from which to stream out data");

  declareProperty("TupleName",
                  m_tupleName = "atlas_ntuple",
                  "Name of the output n-tuple");

  declareProperty("OutputFile",
                  m_outputName,
                  "Name of the output file");
}

RootOutputStreamTool::~RootOutputStreamTool() 
{}

StatusCode 
RootOutputStreamTool::initialize()
{
  ATH_MSG_INFO("Initializing " << name() << " - package version " << PACKAGE_VERSION);

  if (!::AthAlgTool::initialize().isSuccess()) {
    ATH_MSG_FATAL("Cannot initialize AlgTool base class.");
    return(StatusCode::FAILURE);
  }
  // Get the ClassID service
  if (!m_clidSvc.retrieve().isSuccess()) {
    ATH_MSG_FATAL("Cannot get ClassID service via IClassIDSvc interface.");
    return(StatusCode::FAILURE);
  } else {
    ATH_MSG_DEBUG("Found ClassID service.");
  }
  // Get the conversion service
  if (!m_conversionSvc.retrieve().isSuccess()) {
    ATH_MSG_FATAL("Cannot get conversion service via IConversionSvc interface.");
    return(StatusCode::FAILURE);
  } else {
    ATH_MSG_DEBUG("Found conversion service.");
  }
  return(StatusCode::SUCCESS);
}

StatusCode
RootOutputStreamTool::finalize()
{
  // Release the data store service
  if (m_storeSvc != 0) {
    if (!m_storeSvc.release().isSuccess()) {
      ATH_MSG_WARNING("Could not release " << m_storeSvc.type() << " store.");
    }
  }
  // Release the conversion service
  if (!m_conversionSvc.release().isSuccess()) {
    ATH_MSG_WARNING("Cannot release conversion service.");
  }
  // Release the ClassID service
  if (!m_clidSvc.release().isSuccess()) {
    ATH_MSG_WARNING("Cannot release ClassID service.");
  }
  return(::AthAlgTool::finalize());
}

StatusCode
RootOutputStreamTool::connectServices(const std::string& dataStore,
                                      const std::string& cnvSvc,
                                      bool extendProvenenceRecord) 
{
  ATH_MSG_VERBOSE("RootOutputStreamTool::connectServices dataStore = " 
                  << dataStore << ", cnvSvc = " << cnvSvc
                  << ", extendProv=" << extendProvenenceRecord);
  // Release the old data store service
  if (m_storeSvc != 0) {
    if (!m_storeSvc.release().isSuccess()) {
      ATH_MSG_WARNING("Could not release " << m_storeSvc.type() << " store.");
    }
  }
  m_storeSvc = ServiceHandle<StoreGateSvc>(dataStore, this->name());
  // Get the data store service
  if (!m_storeSvc.retrieve().isSuccess()) {
    ATH_MSG_ERROR("Cannot get data store service.");
    return(StatusCode::FAILURE);
  } else {
    ATH_MSG_DEBUG("Found data store service.");
  }
  return(StatusCode::SUCCESS);
}

StatusCode
RootOutputStreamTool::connectOutput(const std::string& outputName) 
{
  ATH_MSG_VERBOSE("RootOutputStreamTool::connectOutput outputName = [" 
                  << outputName <<"]");
  // Set output file name property
  if (!outputName.empty()) {
    m_outputName = outputName;
  } else {
    return(StatusCode::FAILURE);
  }
  // open the file thru the i/o svc
  ServiceHandle<IIoSvc> iosvc("IoSvc/AthIoSvc", name());
  if (!iosvc.retrieve().isSuccess()) {
    ATH_MSG_ERROR("could not retrieve the AthIoSvc");
    return StatusCode::FAILURE;
  }
  IIoSvc::Fd fd = iosvc->open(outputName, IIoSvc::RECREATE);
  if (fd < 0) {
    ATH_MSG_ERROR("could not open-recreate file [" << outputName << "]");
    return StatusCode::FAILURE;
  }
  // FIXME: a better mechanism should be devised to push the tuple-name
  // down to the Athena::RootConnection!!!
  ServiceHandle<IRootSvc> rootsvc("Athena::RootSvc/AthenaRootSvc", name());
  if (!rootsvc.retrieve().isSuccess()) {
    ATH_MSG_ERROR("could not retrieve the AthenaRootSvc");
    return StatusCode::FAILURE;
  }
  // create a RootConnection
  {
    Athena::RootSvc* rsvc = dynamic_cast<Athena::RootSvc*>(&*rootsvc);
    if (!rsvc) {
      ATH_MSG_ERROR("could not dyn-cast to Athena::RootSvc");
      return StatusCode::FAILURE;
    }
    Athena::RootConnection* conn = rsvc->new_connection(fd);
    if (conn == NULL) {
      ATH_MSG_ERROR("could not create a new connection for fd=" << fd 
		    << " fname=[" << outputName << "]");
      return StatusCode::FAILURE;
    }
    conn->setTreeName(m_tupleName);
  }
  // Connect the output file to the service
  if (!m_conversionSvc->connectOutput(m_outputName).isSuccess()) {
    ATH_MSG_ERROR("Unable to connect output " << m_outputName);
    return(StatusCode::FAILURE);
  } else {
    ATH_MSG_DEBUG("Connected to " << m_outputName);
  }
  return(StatusCode::SUCCESS);
}

StatusCode
RootOutputStreamTool::commitOutput() 
{
  ATH_MSG_VERBOSE("RootOutputStreamTool::commitOutput");
  if (m_outputName.empty()) {
    ATH_MSG_ERROR("Unable to commit, no output connected.");
    return(StatusCode::FAILURE);
  }
  // Connect the output file to the service
  if (!m_conversionSvc->commitOutput(m_outputName, false).isSuccess()) {
    ATH_MSG_ERROR("Unable to commit output " << m_outputName);
    return(StatusCode::FAILURE);
  } else {
    ATH_MSG_DEBUG("Committed: " << m_outputName);
  }
  m_outputName.clear();
  return(StatusCode::SUCCESS);
}

StatusCode
RootOutputStreamTool::finalizeOutput() 
{
  ATH_MSG_VERBOSE("RootOutputStreamTool::finalizeOutput");
  return(StatusCode::SUCCESS);
}

StatusCode
RootOutputStreamTool::streamObjects(const IAthenaOutputStreamTool::TypeKeyPairs& typeKeys) 
{
  ATH_MSG_VERBOSE("RootOutputStreamTool::streamObjects(type/keys)...");
  // Now iterate over the type/key pairs and stream out each object
  std::vector<DataObject*> dataObjects;
  dataObjects.reserve(typeKeys.size());
  for (IAthenaOutputStreamTool::TypeKeyPairs::const_iterator 
         first = typeKeys.begin(), 
         last = typeKeys.end();
	   first != last; 
       ++first) {
    const std::string& type = (*first).first;
    const std::string& key  = (*first).second;
    // Find the clid for type name from the classIDSvc
    CLID clid = 0;
    if (!m_clidSvc->getIDOfTypeName(type, clid).isSuccess()) {
      ATH_MSG_ERROR("Could not get clid for typeName " << type);
      return(StatusCode::FAILURE);
    }
    DataObject* dObj = 0;
    // Two options: no key or explicit key
    if (key.empty()) {
      ATH_MSG_DEBUG("Get data object with no key");
      // Get DataObject without key
      dObj = m_storeSvc->accessData(clid);
    } else {
      ATH_MSG_DEBUG("Get data object with key");
      // Get DataObjects with key
      dObj = m_storeSvc->accessData(clid, key);
    }
    if (dObj == 0) {
      // No object - print warning and continue with next object
      ATH_MSG_WARNING("No object found for type " << type << " key " << key);
      continue;
    } else {
      ATH_MSG_DEBUG("Found object for type " << type << " key " << key);
    }
    // Save the dObj
    dataObjects.push_back(dObj);
  }
  return(this->streamObjects(dataObjects));
}

StatusCode
RootOutputStreamTool::streamObjects(const DataObjectVec& dataObjects) 
{
  ATH_MSG_VERBOSE("RootOutputStreamTool::streamObjects(dobjs)");
  if (m_outputName.empty()) {
    ATH_MSG_ERROR("Unable to commit, no output connected.");
    return(StatusCode::FAILURE);
  }
  ATH_MSG_VERBOSE("streaming out... [" << m_conversionSvc.typeAndName() << "]");
  SG::unordered_set<DataObject*> written;
  for (std::vector<DataObject*>::const_iterator 
         doIter = dataObjects.begin(), 
         doLast = dataObjects.end();
       doIter != doLast; 
       ++doIter) {
    ATH_MSG_VERBOSE(" --> [" << (*doIter)->clID() << "/" 
		    << (*doIter)->name() << "]...");
    // Do not stream out same object twice
    if (written.find(*doIter) != written.end()) {
      // Print warning and skip
      ATH_MSG_DEBUG("Trying to write DataObject twice (clid/key): " 
                    << (*doIter)->clID() << ", " << (*doIter)->name());
      ATH_MSG_DEBUG("    Skipping this one.");
    } else {
      written.insert(*doIter);
      // Write object
      IOpaqueAddress* addr = NULL;
      ATH_MSG_VERBOSE(" ==cnvsvc::createRep==...");
      if ((m_conversionSvc->createRep(*doIter, addr)).isSuccess()) {
        IRegistry *ireg = (*doIter)->registry();
        // FIXME: that's a wile hack to handle RootBranchAddress's stickyness.
        // if the transient address already has a RootBranchAddress (ie: it was
        // read from a d3pd-file) calling setAddress will invalidate the
        // previous RBA which will scree things up the next time around.
        // The real fix would probably involve changing a bit the logic either
        // in convsvc::createRep above, or -rather- have the side effect of
        // properly setting up the RootConnection (in createRep) be explicitly
        // written somewhere. (here?)
        if (dynamic_cast<RootBranchAddress*>(ireg->address())) {
          delete addr; addr = NULL;
        } else {
          ireg->setAddress(addr);
        }
        // SG::DataProxy* proxy = dynamic_cast<SG::DataProxy*>((*doIter)->registry());
        // if (!proxy) {
        //   ATH_MSG_WARNING("Could cast DataObject " 
        //                   << (*doIter)->clID() << " " << (*doIter)->name());
        // }
      } else {
        ATH_MSG_ERROR("Could not create Rep for DataObject (clid/key): " 
                      << (*doIter)->clID() << ", " << (*doIter)->name());
        return(StatusCode::FAILURE);
      }
    }
  }
  return(StatusCode::SUCCESS);
}

StatusCode
RootOutputStreamTool::fillObjectRefs(const DataObjectVec& /*dataObjects*/) 
{
  ATH_MSG_VERBOSE("RootOutputStreamTool::fillObjectRefs");
  return(StatusCode::SUCCESS);
}

StatusCode 
RootOutputStreamTool::getInputItemList(SG::IFolder* /*m_p2BWrittenFromTool*/) 
{
  ATH_MSG_VERBOSE("RootOutputStreamTool::getInputItemList");
  return(StatusCode::SUCCESS);
}

}//> namespace Athena
