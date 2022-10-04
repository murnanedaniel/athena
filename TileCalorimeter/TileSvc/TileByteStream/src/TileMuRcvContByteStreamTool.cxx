/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//****************************************************************************
// Filename : TileMuRcvContByteStreamTool.cxx
// Author   : Joao Gentil Saraiva (jmendes@cern.ch)
// Created  : February 2015
//
// DESCRIPTION
//    AlgTool used in the TileMuRcv to BS conversion
//
// BUGS:
//
// History:
//
//****************************************************************************

// Gaudi includes

// Atlas includes
#include "AthenaKernel/errorcheck.h"
#include "StoreGate/ReadCondHandle.h"

// Tile includes
#include "TileByteStream/TileMuRcvContByteStreamTool.h"
#include "TileByteStream/TileROD_Decoder.h"
#include "TileByteStream/TileROD_Encoder.h"
#include "TileEvent/TileMuonReceiverObj.h"
#include "TileEvent/TileContainer.h"
#include "TileIdentifier/TileHWID.h"

#include "AthenaKernel/CLASS_DEF.h"

#include <map> 
#include <stdint.h>

static const InterfaceID IID_ITileMuRcvContByteStreamTool("TileMuRcvContByteStreamTool", 1, 0);

const InterfaceID& TileMuRcvContByteStreamTool::interfaceID() {
  return IID_ITileMuRcvContByteStreamTool;
}

// default constructor

TileMuRcvContByteStreamTool::TileMuRcvContByteStreamTool(const std::string& type, const std::string& name,
    const IInterface* parent)
  : AthAlgTool(type, name, parent)
  , m_tileHWID(0)
  , m_runPeriod(0)
{
  declareInterface<TileMuRcvContByteStreamTool>(this);
}

// destructor

TileMuRcvContByteStreamTool::~TileMuRcvContByteStreamTool() {
}

StatusCode TileMuRcvContByteStreamTool::initialize() {

  ATH_MSG_INFO ("Initializing TileMuRcvContByteStreamTool");

  ATH_CHECK( detStore()->retrieve(m_tileHWID, "TileHWID") );

  ToolHandle<TileROD_Decoder> dec("TileROD_Decoder");
  ATH_CHECK( dec.retrieve() );

  ATH_CHECK( m_hid2RESrcIDKey.initialize(m_initializeForWriting) );

  const TileCablingService *cabling = TileCablingService::getInstance();
  m_runPeriod = cabling->runPeriod();

  return StatusCode::SUCCESS;
}

StatusCode TileMuRcvContByteStreamTool::finalize() {
  ATH_MSG_INFO ("Finalizing TileMuRcvContByteStreamTool successfuly");
  return StatusCode::SUCCESS;
}

StatusCode TileMuRcvContByteStreamTool::convert(TileMuonReceiverContainer* cont, FullEventAssembler<TileHid2RESrcID> *fea) const
{
  ATH_MSG_INFO ("Executing TileMuRcvContByteStreamTool::convert method");

  SG::ReadCondHandle<TileHid2RESrcID> hid2re{m_hid2RESrcIDKey};

  int  n           = 0;
  uint32_t frag_id = 0x0;
  uint32_t reid    = 0x0;
  
  TileMuonReceiverContainer::const_iterator it_cont  = cont->begin();
  TileMuonReceiverContainer::const_iterator end_cont = cont->end();

  // skip thresholds stored at first position of the container
  //
  ++it_cont;

  std::map<uint32_t, TileROD_Encoder> mapEncoder;

  for (; it_cont != end_cont; ++it_cont) 
    {
      n++;
      frag_id = (*it_cont)->identify();
      reid = hid2re->getRodTileMuRcvID(frag_id);
      mapEncoder[reid].setTileHWID(m_tileHWID,m_runPeriod);
      const TileMuonReceiverObj* tileMuRcv = *it_cont;	
      mapEncoder[reid].addTileMuRcvObj(tileMuRcv);
    }                                                            

  ATH_MSG_DEBUG( " Number of TileMuonReceiverObj objects counted " << n << " out of the possible " << cont->size()-1 ); 

  // fill ROD
  //
  std::map<uint32_t, TileROD_Encoder>::iterator map_it  = mapEncoder.begin();
  std::map<uint32_t, TileROD_Encoder>::iterator map_end = mapEncoder.end();

  FullEventAssembler<TileHid2RESrcID>::RODDATA* theROD;

  TileROD_Encoder* theEncoder;

  for (; map_it != map_end; ++map_it) 
    {
      theROD     = fea->getRodData( (*map_it).first );
      theEncoder = &( (*map_it).second );
      theEncoder -> fillRODTileMuRcvObj( *theROD );
      ATH_MSG_DEBUG( " Number of words in ROD " <<MSG::hex<< (*map_it).first <<MSG::dec<< ": " << theROD->size() );// DEBUG
      ATH_MSG_DEBUG( " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " );// DEBUG
    }

  return StatusCode::SUCCESS;
}
