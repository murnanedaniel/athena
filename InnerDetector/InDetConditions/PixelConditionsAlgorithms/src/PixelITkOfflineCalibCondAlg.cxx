/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#include "PixelITkOfflineCalibCondAlg.h"
#include "Identifier/Identifier.h"
#include "GaudiKernel/EventIDRange.h"
#include "PathResolver/PathResolver.h"
#include <memory>
#include <sstream>


PixelITkOfflineCalibCondAlg::PixelITkOfflineCalibCondAlg(const std::string& name, ISvcLocator* pSvcLocator):
  ::AthReentrantAlgorithm(name, pSvcLocator),
  m_condSvc("CondSvc", name),
  // For 21.9 only, to be removed in master
  m_detStore(nullptr)
{
  declareProperty("InputSource",m_inputSource=2,"Source of data: 0 (none), 1 (text file), 2 (database)");
  declareProperty("PixelClusterErrorDataFile", m_textFileName="PixelITkClusterErrorData.txt","Read constants from this file");  
  declareProperty("DumpConstants", m_dump=0, "Dump constants to text file"); 
}



StatusCode PixelITkOfflineCalibCondAlg::initialize() {
  ATH_MSG_DEBUG("PixelITkOfflineCalibCondAlg::initialize()");

  ATH_CHECK(m_condSvc.retrieve());

  // For 21.9 only, to be removed in master
  ATH_CHECK(service("DetectorStore", m_detStore));

  if (m_inputSource==2){
    if(m_readKey.key().empty()) {
     ATH_MSG_ERROR("The database is set to be input source (2) but the ReadKey is empty.");
    }
    ATH_CHECK(m_readKey.initialize());
  }

  ATH_CHECK(m_writeKey.initialize());

  if (m_condSvc->regHandle(this,m_writeKey,m_writeKey.key()).isFailure()) {
    ATH_MSG_FATAL("unable to register WriteCondHandle " << m_writeKey.fullKey() << " with CondSvc");
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}





StatusCode PixelITkOfflineCalibCondAlg::execute_r(const EventContext& ctx) const {
  ATH_MSG_DEBUG("PixelITkOfflineCalibCondAlg::execute_r()");

  SG::WriteCondHandle<PixelCalib::PixelITkOfflineCalibData> writeHandle(m_writeKey, ctx);
  if (writeHandle.isValid()) {
    ATH_MSG_DEBUG("CondHandle " << writeHandle.fullKey() << " is already valid.. In theory this should not be called, but may happen if multiple concurrent events are being processed out of order.");
    return StatusCode::SUCCESS; 
  }

  // For 21.9 only, to be removed in master
  bool inDetectorStore = m_detStore->contains<PixelCalib::PixelITkOfflineCalibData>(writeHandle.key());
  if(inDetectorStore) {
    ATH_MSG_DEBUG(writeHandle.fullKey() << " already in DetectorStore");
    return StatusCode::SUCCESS;
  }


  // Construct the output Cond Object and fill it in
  PixelCalib::PixelITkOfflineCalibData* writeCdo = new PixelCalib::PixelITkOfflineCalibData();

  if (m_inputSource==0) {
    ATH_MSG_WARNING("So far do nothing!! return StatusCode::FAILURE");
    return StatusCode::FAILURE;
  }
  else if (m_inputSource==1) {
    ATH_MSG_WARNING("Pixel ITk constants read from text file. Only supported for local developments and debugging!");

    PixelCalib::PixelITkOfflineCalibData* calibData = new PixelCalib::PixelITkOfflineCalibData;
    
    PixelCalib::PixelITkClusterErrorData* pced = calibData->getPixelITkClusterErrorData();

    // Find and open the text file
    ATH_MSG_INFO("Load PixelITkErrorData constants from text file");
    std::string fileName = PathResolver::find_file(m_textFileName, "DATAPATH");
    if (fileName.size()==0) { ATH_MSG_WARNING("Input file " << fileName << " not found! Default (hardwired) values to be used!"); }
    else { pced->load(fileName);  }

    ATH_MSG_DEBUG("Get error constants");
    std::vector<float> constants = calibData->getConstants();
    if (constants.size()) { ATH_MSG_VERBOSE("constants are defined"); }
    else                  { ATH_MSG_ERROR("constants size is NULL!!!"); } 


    const EventIDBase start{EventIDBase::UNDEFNUM, EventIDBase::UNDEFEVT, 0,                       0,                       EventIDBase::UNDEFNUM, EventIDBase::UNDEFNUM};
    const EventIDBase stop {EventIDBase::UNDEFNUM, EventIDBase::UNDEFEVT, EventIDBase::UNDEFNUM-1, EventIDBase::UNDEFNUM-1, EventIDBase::UNDEFNUM, EventIDBase::UNDEFNUM};
    const EventIDRange rangeW{start, stop};

    ATH_MSG_DEBUG("Range of input is " << rangeW);
    
    if (constants.size()) {
      ATH_MSG_DEBUG("Found constants with new-style Identifier key");
      writeCdo->setConstants(constants);
    }

    if (writeHandle.record(rangeW, writeCdo).isFailure()) {
      ATH_MSG_FATAL("Could not record PixelCalib::PixelITkOfflineCalibData " << writeHandle.key() << " with EventRange " << rangeW << " into Conditions Store");
      return StatusCode::FAILURE;
    }

    ATH_MSG_DEBUG("recorded new CDO " << writeHandle.key() << " with range " << rangeW << " into Conditions Store");


    // For 21.9 only, to be removed in master

    if(m_detStore->record(writeCdo,writeHandle.key()).isFailure()) {
      ATH_MSG_FATAL("Could not record PixelCalib::PixelITkOfflineCalibData " << writeHandle.key() << "into StoreGate");
      return StatusCode::FAILURE;
    }

    if (m_dump!=0) {
      ATH_MSG_DEBUG("Dump the constants to file");
      calibData->dump();
    }
    delete calibData;

  }

  else if (m_inputSource==2) {

    SG::ReadCondHandle<CondAttrListCollection> readHandle{m_readKey, ctx};

    const EventIDBase start{EventIDBase::UNDEFNUM, EventIDBase::UNDEFEVT, 0,                       0,                       EventIDBase::UNDEFNUM, EventIDBase::UNDEFNUM};
    const EventIDBase stop {EventIDBase::UNDEFNUM, EventIDBase::UNDEFEVT, EventIDBase::UNDEFNUM-1, EventIDBase::UNDEFNUM-1, EventIDBase::UNDEFNUM, EventIDBase::UNDEFNUM};
    EventIDRange rangeW{start, stop};

    // To be used in master

    /*const CondAttrListCollection* readCdo{*readHandle};

    // Get the validitiy range
    if (not readHandle.range(rangeW)) {
       ATH_MSG_FATAL("Failed to retrieve validity range for " << readHandle.key());
       return StatusCode::FAILURE;
    }*/

    // Temporary workaround in 21.9

    const CondAttrListCollection* readCdo;
    if(m_detStore->retrieve(readCdo,readHandle.key()).isFailure()) {
      ATH_MSG_FATAL("Could not retrieve CondAttrListCollection " << readHandle.key() << " from StoreGate");
      return StatusCode::FAILURE;
    }

    if (readCdo==nullptr) {
      ATH_MSG_FATAL("Null pointer to the read conditions object");
      return StatusCode::FAILURE;
    }

    ATH_MSG_DEBUG("Range of input is " << rangeW);

    std::vector<float> constants;

    for(CondAttrListCollection::const_iterator attrList = readCdo->begin(); attrList != readCdo->end(); ++attrList){

      std::ostringstream attrStr;
      (*attrList).second.toOutputStream(attrStr);
      ATH_MSG_DEBUG( "ChanNum " << (*attrList).first << " Attribute list " << attrStr.str() );

      constants.emplace_back( (*attrList).second["pixelID"].data<long long>() );
      constants.emplace_back( (*attrList).second["delta_x"].data<float>() );
      constants.emplace_back( (*attrList).second["delta_error_x"].data<float>() );
      constants.emplace_back( (*attrList).second["delta_y"].data<float>() );
      constants.emplace_back( (*attrList).second["delta_error_y"].data<float>() );

    }

    writeCdo->setConstants(constants);


    if (writeHandle.record(rangeW, writeCdo).isFailure()) {
      ATH_MSG_FATAL("Could not record PixelCalib::PixelITkOfflineCalibData " << writeHandle.key() << " with EventRange " << rangeW << " into Conditions Store");
      return StatusCode::FAILURE;
    }
    ATH_MSG_DEBUG("recorded new CDO " << writeHandle.key() << " with range " << rangeW << " into Conditions Store");

    // For 21.9 only, to be removed in master

    if(m_detStore->record(writeCdo,writeHandle.key()).isFailure()) {
      ATH_MSG_FATAL("Could not record PixelCalib::PixelITkOfflineCalibData " << writeHandle.key() << "into StoreGate");
      return StatusCode::FAILURE;
    }

  }

  return StatusCode::SUCCESS;
}


StatusCode PixelITkOfflineCalibCondAlg::finalize() {
  ATH_MSG_DEBUG("PixelITkOfflineCalibCondAlg::finalize()");
  return StatusCode::SUCCESS;
}

