/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "LArBadChannelTool/LArBadFebCondAlg.h"
#include "LArBadChannelTool/LArBadChanBlobUtils.h"
#include "LArIdentifier/LArOnlineID.h"
#include "LArBadChannelTool/LArBadChannelDecoder.h"
#include "AthenaKernel/IOVInfiniteRange.h"

StatusCode LArBadFebCondAlg::initialize() {

  // Read Handles
  ATH_CHECK( m_BCInputKey.initialize(SG::AllowEmpty) );
  ATH_CHECK( m_BCOutputKey.initialize() );

  return StatusCode::SUCCESS;
}


StatusCode LArBadFebCondAlg::execute(const EventContext& ctx) const {
    
  SG::WriteCondHandle<LArBadFebCont> writeHandle{m_BCOutputKey,ctx};
  
  if (writeHandle.isValid()) {
    msg(MSG::DEBUG) << "Found valid write handle" << endmsg;
    return StatusCode::SUCCESS;
  }  

  std::unique_ptr<LArBadFebCont> badFebCont=std::make_unique<LArBadFebCont>();


  if(!m_BCInputKey.key().empty()) {

    SG::ReadCondHandle<AthenaAttributeList> readHandle{m_BCInputKey,ctx};
    const AthenaAttributeList* attrList{*readHandle};
    if (attrList==nullptr) {
      msg(MSG::ERROR) << "Failed to retrieve CondAttributeListCollection with key " << m_BCInputKey.key() << endmsg;
      return StatusCode::FAILURE;
    }
    writeHandle.addDependency(readHandle); 
    
    const coral::Blob& blob = (*attrList)["Blob"].data<coral::Blob>();
    unsigned int chanSize = (*attrList)["ChannelSize"].data<unsigned int>();
    unsigned int stateSize = (*attrList)["StatusWordSize"].data<unsigned int>();
    unsigned int endian = (*attrList)["Endianness"].data<unsigned int>();
    unsigned int version = (*attrList)["Version"].data<unsigned int>();
      
    std::vector<std::pair<HWIdentifier,LArBadFeb> > bcVec = 
      LArBadChanBlobUtils::decodeBlob<LArBadFeb>( &blob, chanSize, stateSize, endian,
          					version, msg());
     
    for (auto& idBC : bcVec) {
      badFebCont->add(idBC.first,idBC.second);
    }
  }
  else {
    //No input data from DB, set infinte range
    writeHandle.addDependency(IOVInfiniteRange::infiniteRunLB());
  }

  if (m_inputFileName.size()) {//Read supplemental data from ASCII file (if required)
     
     const LArOnlineID* onlineID;
     ATH_CHECK(detStore()->retrieve(onlineID,"LArOnlineID"));	       
     LArBadChannelDecoder decoder(&(*onlineID));
     std::vector<std::pair<HWIdentifier,LArBadFeb> > bcVec = decoder.readFebASCII(m_inputFileName, msg());
     for (auto& idBC : bcVec) {
       badFebCont->add(idBC.first,idBC.second);
     }
  } //end if have ASCII filename
 
 
 
  size_t nChanBeforeMege=badFebCont->size();
  badFebCont->sort(); //Sorts vector of bad febs and merges duplicate entries
   
  ATH_MSG_INFO("Read a total of " << badFebCont->size() << " problematic febs from database");
  if (nChanBeforeMege!=badFebCont->size()) {
    ATH_MSG_INFO("Merged " << nChanBeforeMege-badFebCont->size() << " duplicate entries");
  }
 
  if(writeHandle.record(std::move(badFebCont)).isFailure()) {
    ATH_MSG_ERROR("Could not record LArBadFebCont object with " 
		  << writeHandle.key() 
		  << " with EventRange " << writeHandle.getRange()
		  << " into Conditions Store");
    return StatusCode::FAILURE;
  }
  ATH_MSG_INFO("Recorded LArBadFebCont object with " 
	       << writeHandle.key() 
	       << " with EventRange " << writeHandle.getRange()
	       << " into Conditions Store");


  return StatusCode::SUCCESS;
}
