/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


// Gaudi/Athena include(s):
#include "GaudiKernel/MsgStream.h"

// CTP data format include(s):
#include "CTPfragment/CTPfragment.h"
#include "CTPfragment/CTPdataformat.h"
#include "L1CommonCore/CTPdataformatVersion.h"

// Trigger include(s):
#include "TrigT1Result/CTP_RDO.h"

// Local include(s):
#include "TrigT1ResultByteStream/CTPByteStreamTool.h"

/// Unique interface ID of the tool that identifies it to the framweork
static const InterfaceID IID_ICTPByteStreamTool( "CTPByteStreamTool", 1, 0 );

/**
 * This function is needed by the framework to know what kind of tool
 * this is.
 */
const InterfaceID& CTPByteStreamTool::interfaceID() {
  return IID_ICTPByteStreamTool;
}

/**
 * The constructor takes care of correctly constructing the base class and
 * declaring the tool's interface to the framework.
 */
CTPByteStreamTool::CTPByteStreamTool( const std::string& type, const std::string& name,
                                      const IInterface* parent )
  : AlgTool( type, name, parent ), m_srcIdMap( 0 ) {
  
  declareInterface< CTPByteStreamTool >( this );
}

/**
 * The destructor doesn't do anything.
 */
CTPByteStreamTool::~CTPByteStreamTool() {

}

/**
 * The function creates a CTPSrcIdMap object that is used in the conversion
 * and initialises the base class.
 */
StatusCode CTPByteStreamTool::initialize() {

  m_srcIdMap = new CTPSrcIdMap();
  return AlgTool::initialize();

}

/**
 * The function deletes the CTPSrcIdMap object and finalises the base class.
 */
StatusCode CTPByteStreamTool::finalize() {

  delete m_srcIdMap;
  return AlgTool::finalize();

}


/**
 * Conversion from RDO to eformat::ROBFragment.
 * This is called from the CTPByteStreamCnv::createRep method.
 */
StatusCode CTPByteStreamTool::convert( const CTP_RDO* result, RawEventWrite* re ) {

  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "executing convert() from RDO to ROBFragment" << endreq;

  //get CTP version
  unsigned int ctpVersionNumber = result->getCTPVersionNumber();
  CTPdataformatVersion ctpVersion(ctpVersionNumber);
  
  // Clear Event Assembler
  m_fea.clear();

  // Reset lumi-block number to 1
  m_fea.setDetEvtType( 1 );
  // Set L1Apos to center of readout window
  uint16_t minorVersion = ( result->getNumberOfBunches() - 1u ) / 2u;
  minorVersion &= ctpVersion.getL1APositionMask();
  minorVersion <<= ctpVersion.getL1APositionShift();
  m_fea.setRodMinorVersion( minorVersion );  

  FullEventAssembler< CTPSrcIdMap >::RODDATA* theROD;

  // Source ID of CTP
  const uint32_t rodId = m_srcIdMap->getRodID();

  log << MSG::DEBUG << " ROD ID:" << MSG::hex << rodId << endreq;

  // get the ROD data container to be filled
  theROD = m_fea.getRodData( rodId );

  log << MSG::VERBOSE << " Dumping CTP words:" << endreq;

  // fill Data Words
  const std::vector< uint32_t >& vDataWords = result->getDataWords();
  std::vector< uint32_t >::const_iterator it = vDataWords.begin();
  std::vector< uint32_t >::const_iterator end = vDataWords.end();
  for( ; it != end; ++it ) {
    theROD->push_back( *it );
    log << MSG::VERBOSE << "     0x" << MSG::hex << std::setw( 8 )  << ( *it ) << endreq;
  }

  // Now fill full event
  log << MSG::DEBUG << "Now filling the event with the CTP fragment" << endreq;
  m_fea.fill( re, log );

  return StatusCode::SUCCESS;
}


/**
 * Conversion from eformat::ROBFragment to RDO.
 * This is called from the CTPByteStreamCnv::createObj method.
 */
StatusCode CTPByteStreamTool::convert( const ROBF* rob, CTP_RDO*& result ) {

  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "executing convert() from ROBFragment to RDO" << endreq;
  
  const uint32_t ctpRodId = m_srcIdMap->getRodID();
  const uint32_t rodId = rob->rod_source_id();

  log << MSG::DEBUG << " expected ROD sub-detector ID: " << std::hex << ctpRodId 
      << " ID found: " << std::hex << rodId << std::dec << endreq;  

  if ( rodId == ctpRodId ) {

    log << MSG::VERBOSE << " Dumping CTP words:" << endreq;

    // For generality let's declare the data pointer like this. Altough it's
    // unlikely to ever change from being a pointer to uint32_t-s.
    OFFLINE_FRAGMENTS_NAMESPACE::PointerType it_data;
    rob->rod_data( it_data );
    const uint32_t ndata = rob->rod_ndata();

    // data words
    std::vector<uint32_t> vDataWords;
    for ( uint32_t i = 0; i < ndata; ++i, ++it_data ) {
      vDataWords.push_back( static_cast<uint32_t>( *it_data ) );
      log << MSG::VERBOSE << "     0x" << MSG::hex << std::setw( 8 )  << *it_data << endreq;
    }

    // create CTP RDO
    uint32_t nExtraWords=0;
    nExtraWords=CTPfragment::numberExtraPayloadWords(rob);
    unsigned int ctpVersionNumber = CTPfragment::ctpFormatVersion(rob);
    result = new CTP_RDO(ctpVersionNumber, vDataWords , nExtraWords);

    uint8_t l1apos =  CTPfragment::lvl1AcceptBunch(rob);
    result->setL1AcceptBunchPosition(l1apos);

    uint32_t tc = CTPfragment::turnCounter(rob);
    result->setTurnCounter(tc);

    return StatusCode::SUCCESS;
  }

  log << MSG::ERROR << "Wrong ROD ID found in the CTP ROB fragment!" << endreq;
  return StatusCode::FAILURE;
}
