/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/


#include "AthenaKernel/errorcheck.h"
#include "StoreGate/ReadCondHandle.h"

#include "TileByteStream/TileRawChannelContByteStreamTool.h"
#include "TileByteStream/TileROD_Encoder.h"

#include "TileEvent/TileRawChannelCollection.h"
#include "TileEvent/TileRawChannelContainer.h"
#include "TileEvent/TileRawChannel.h"
#include "TileEvent/TileFastRawChannel.h"

#include "TileIdentifier/TileHWID.h"
#include "TileCalibBlobObjs/TileCalibUtils.h"
#include "TileConditions/TileCondToolEmscale.h"
#include "TileConditions/ITileBadChanTool.h"
#include "TileConditions/TileCablingService.h"
#include "TileByteStream/TileROD_Decoder.h"

#include <map> 
#include <stdint.h>

static const InterfaceID IID_ITileRawChannelContByteStreamTool("TileRawChannelContByteStreamTool", 1, 0);

const InterfaceID& TileRawChannelContByteStreamTool::interfaceID() {
  return IID_ITileRawChannelContByteStreamTool;
}

// default constructor

TileRawChannelContByteStreamTool::TileRawChannelContByteStreamTool(const std::string& type,
    const std::string& name, const IInterface* parent)
    : AthAlgTool(type, name, parent)
    , m_tileHWID(0)
    , m_verbose(false)
    , m_maxChannels(TileCalibUtils::MAX_CHAN)
{
  declareInterface<TileRawChannelContByteStreamTool>(this);
}

// destructor 

TileRawChannelContByteStreamTool::~TileRawChannelContByteStreamTool() {
}

StatusCode TileRawChannelContByteStreamTool::initialize() {

  ATH_MSG_INFO ("Initializing TileRawChannelContByteStreamTool");

  ATH_CHECK( detStore()->retrieve(m_tileHWID, "TileHWID") );

  ToolHandle<TileROD_Decoder> dec("TileROD_Decoder");
  ATH_CHECK( dec.retrieve() );

  // get TileCondToolEmscale
  ATH_CHECK( m_tileToolEmscale.retrieve() );

  // get TileBadChanTool
  ATH_CHECK( m_tileBadChanTool.retrieve() );

  ATH_CHECK( m_hid2RESrcIDKey.initialize(m_initializeForWriting) );

  m_maxChannels = TileCablingService::getInstance()->getMaxChannels();

  return StatusCode::SUCCESS;
}

StatusCode TileRawChannelContByteStreamTool::finalize() {
  ATH_MSG_INFO ("Finalizing TileRawChannelContByteStreamTool successfuly");
  return StatusCode::SUCCESS;
}

StatusCode TileRawChannelContByteStreamTool::convert(CONTAINER* rawChannelContainer, FullEventAssembler<TileHid2RESrcID> *fea) const
{
  bool isTMDB = evtStore()->proxy(rawChannelContainer)->name() == "MuRcvRawChCnt";

  TileFragHash::TYPE contType = rawChannelContainer->get_type();
  TileRawChannelUnit::UNIT inputUnit = rawChannelContainer->get_unit();
  TileRawChannelUnit::UNIT outputUnit = inputUnit; 

  bool oflCont = (inputUnit < TileRawChannelUnit::OnlineOffset);

  FullEventAssembler<TileHid2RESrcID>::RODDATA* theROD;
  SG::ReadCondHandle<TileHid2RESrcID> hid2re{m_hid2RESrcIDKey};

  ATH_MSG_DEBUG( " Number of raw channel collections... " << rawChannelContainer->size() << " " << evtStore()->proxy(rawChannelContainer)->name());

  std::map<uint32_t, TileROD_Encoder> mapEncoder;
  std::vector<TileFastRawChannel> channels;
  channels.reserve (m_tileHWID->channel_hash_max());

  int nch = 0;
  uint32_t reid = 0x0;

  for (const TileRawChannelCollection* rawChannelCollection : *rawChannelContainer) {

    TileRawChannelCollection::ID frag_id = rawChannelCollection->identify();

    if (isTMDB) reid = hid2re->getRodTileMuRcvID(frag_id);
    else reid = hid2re->getRodID(frag_id);

    TileROD_Encoder& encoder = mapEncoder[reid];

    encoder.setTileHWID(m_tileHWID, m_verbose, 4);
    encoder.setTypeAndUnit(contType, outputUnit);
    encoder.setMaxChannels(m_maxChannels);

    HWIdentifier drawer_id = m_tileHWID->drawer_id(frag_id);

    int ros = m_tileHWID->ros(drawer_id);
    int drawer = m_tileHWID->drawer(drawer_id);
    int drawerIdx = TileCalibUtils::getDrawerIdx(ros, drawer);

    int n = 0;

    for (const TileRawChannel* rawChannel : *rawChannelCollection) {

      HWIdentifier adc_id = rawChannel->adc_HWID();
      int channel = m_tileHWID->channel(adc_id);
      int adc = m_tileHWID->adc(adc_id);
      float amplitude = rawChannel->amplitude();
      float time = rawChannel->time();
      float quality = rawChannel->quality();
      if (isTMDB) {
        channels.emplace_back (frag_id, channel, adc, amplitude, 0., 0.);
      } else {
        if (oflCont) {
          if (quality > 15.0) quality = 15.0;
          if (m_tileBadChanTool->getAdcStatus(drawerIdx, channel, adc).isBad()) quality += 16.;
	}
        //amplitude = m_tileToolEmscale->channelCalib(drawerIdx, channel, adc, amplitude, inputUnit, outputUnit);
        channels.emplace_back (frag_id, channel, adc, amplitude, time, quality);
      }

      // Don't need to worry about these moving due to the reserve() above.
      encoder.add(&channels.back());
      ++nch;
      ++n;
    }

    ATH_MSG_DEBUG( " Collection " << MSG::hex << "0x" << frag_id
                  << " ROD " << "0x" << reid
                  << " number of channels " << MSG::dec << n );
  }

  // TileROD_Encoder has collected all the channels, now can fill the
  // ROD block data.

  for (std::pair<const uint32_t, TileROD_Encoder>& reidAndEncoder: mapEncoder) {

    theROD = fea->getRodData(reidAndEncoder.first);
    TileROD_Encoder& theEncoder = reidAndEncoder.second;

    if ((reidAndEncoder.first & 0xf00)) {
       theEncoder.fillRODTileMuRcvRawChannel(*theROD);
    } else {
      if (m_doFragType4) theEncoder.fillROD4(*theROD);
      if (m_doFragType5) theEncoder.fillROD5(*theROD);
    }
    ATH_MSG_DEBUG( " Number of TileRawChannel words in ROD " << MSG::hex << " 0x" << reidAndEncoder.first << MSG::dec << " : " << theROD->size() );
  }

  return StatusCode::SUCCESS;
}
