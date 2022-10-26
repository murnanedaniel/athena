/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#include "AFP_Raw2DigiTool.h"
#include "AFP_RawEv/AFP_RawContainer.h"
#include "AFP_RawEv/AFP_ROBID.h"
#include "xAODForward/AFPStationID.h"

#include "xAODForward/AFPSiHit.h"
#include "xAODForward/AFPSiHitContainer.h"
#include "xAODForward/AFPSiHitAuxContainer.h"

#include "xAODForward/AFPToFHit.h"
#include "xAODForward/AFPToFHitContainer.h"
#include "xAODForward/AFPToFHitAuxContainer.h"

#include <vector>

AFP_Raw2DigiTool::AFP_Raw2DigiTool(const std::string &type,
				   const std::string &name,
				   const IInterface *parent)
  : base_class(type, name, parent),
    m_totToChargeTransformation ("totToChargeTransformation", "1909 + x*363 + x*x*141")
{}		

AFP_Raw2DigiTool::~AFP_Raw2DigiTool() {}

StatusCode AFP_Raw2DigiTool::initialize()
{
  ATH_MSG_DEBUG("Initializing " << name() << "...");

  ATH_CHECK( m_rawDataContainerName.initialize() );
  ATH_CHECK( m_AFPSiHitsContainerName.initialize() );
  ATH_CHECK( m_AFPHitsContainerNameToF.initialize() );
  
  m_totToChargeTransformation=TF1(m_totToChargeTransfName.toString().c_str(), m_totToChargeTransfExpr.toString().c_str());

  return StatusCode::SUCCESS;
}

StatusCode AFP_Raw2DigiTool::recoAll(const EventContext &ctx) const
{
  ATH_CHECK( recoSiHits(ctx) );
  ATH_CHECK( recoToFHits(ctx) );

  return StatusCode::SUCCESS;
}

StatusCode AFP_Raw2DigiTool::recoSiHits(const EventContext &ctx) const
{
  ATH_MSG_DEBUG("AFP_Raw2DigiTool recoSiHits ");

  // create output containers
  auto siHitContainer = std::make_unique<xAOD::AFPSiHitContainer>();
  auto siHitAuxContainer = std::make_unique<xAOD::AFPSiHitAuxContainer>();
  siHitContainer->setStore(siHitAuxContainer.get());

  // retrieve raw data
  SG::ReadHandle<AFP_RawContainer> container{m_rawDataContainerName, ctx};
  if (!container.isValid()) {
    ATH_MSG_WARNING("AFP_Raw2DigiTool: Could not find raw data container");
    return StatusCode::SUCCESS;
  }
  else
    ATH_MSG_DEBUG("AFP_Raw2DigiTool: Raw data container retrieved");

  for (const AFP_SiRawCollection& collection: container->collectionsSi())
    for (const AFP_SiRawData& data : collection.dataRecords())
      newXAODHitSi (siHitContainer.get(), collection, data);

  SG::WriteHandle<xAOD::AFPSiHitContainer> writeHandle{m_AFPSiHitsContainerName, ctx};
  ATH_CHECK( writeHandle.record(std::move(siHitContainer), std::move(siHitAuxContainer)) );

  return StatusCode::SUCCESS;
}

StatusCode AFP_Raw2DigiTool::recoToFHits(const EventContext &ctx) const
{
  ATH_MSG_DEBUG("AFP_Raw2DigiTool recoToFHits ");

  // create output containers
  auto tofHitContainer = std::make_unique<xAOD::AFPToFHitContainer>();
  auto tofHitAuxContainer = std::make_unique<xAOD::AFPToFHitAuxContainer>();
  tofHitContainer->setStore(tofHitAuxContainer.get());

  // retrieve raw data
  SG::ReadHandle<AFP_RawContainer> container{m_rawDataContainerName, ctx};
  if (!container.isValid()) {
    ATH_MSG_WARNING("AFP_Raw2DigiTool: Could not find raw data container");
    return StatusCode::SUCCESS;
  }
  else
    ATH_MSG_DEBUG("AFP_Raw2DigiTool: Raw data container retrieved");
  
  for (const AFP_ToFRawCollection& collection: container->collectionsToF())
    for (const AFP_ToFRawData& data : collection.dataRecords())
      if (data.hitDiscConfig() == 3 && (data.header() == 2 || (data.header() == 1 && !data.isTrigger()))) 
	newXAODHitToF (tofHitContainer.get(), collection, data, ctx); // is_ToF && (is_HPTDC || (is_picoTDC && !is_trigger_word))

  SG::WriteHandle<xAOD::AFPToFHitContainer> writeHandle{m_AFPHitsContainerNameToF, ctx};
  ATH_CHECK( writeHandle.record(std::move(tofHitContainer), std::move(tofHitAuxContainer)) );

  return StatusCode::SUCCESS;
}


StatusCode AFP_Raw2DigiTool::finalize(){
  ATH_MSG_INFO("Finalizing " << name() << "...");


  return StatusCode::SUCCESS;
}


void AFP_Raw2DigiTool::newXAODHitToF (xAOD::AFPToFHitContainer* tofHitContainer, const AFP_ToFRawCollection& collection, const AFP_ToFRawData& data, const EventContext& ctx) const
{
  xAOD::AFPToFHit* xAODToFHit = new xAOD::AFPToFHit();
  tofHitContainer->push_back(xAODToFHit);

  const uint32_t hitLink = data.link();
  const unsigned int robID = collection.robId();
  xAODToFHit->setTime(data.time()*s_timeConversionFactor);
  xAODToFHit->setPulseLength(data.pulseLength()*s_pulseLengthFactor); 
  xAODToFHit->setHptdcChannel(data.channel());

  // set station ID
  if (robID == AFP_ROBID::sideC) // no AFP_ROBID::sideC_2016 because there was no ToF in 2016
    xAODToFHit->setStationID(xAOD::AFPStationID::farC);
  else if (robID == AFP_ROBID::sideA)
    xAODToFHit->setStationID(xAOD::AFPStationID::farA);
  else
    xAODToFHit->setStationID(-1);

  // set hptdcID
  if (hitLink == 12)
    xAODToFHit->setHptdcID(2);
  else if (hitLink == 13)
    xAODToFHit->setHptdcID(1);
  else
    xAODToFHit->setHptdcID(-1);

  // set barID
  setBarAndTrainID(xAODToFHit, ctx);
}

void AFP_Raw2DigiTool::newXAODHitSi (xAOD::AFPSiHitContainer* siHitContainer, const AFP_SiRawCollection& collection, const AFP_SiRawData& data) const
{
  xAOD::AFPSiHit* xAODSiHit = new xAOD::AFPSiHit();
  siHitContainer->push_back(xAODSiHit);

  const unsigned int hitLink = data.link();
  const unsigned int robID = collection.robId();
  // decoding according to Outlink numbers in: https://twiki.cern.ch/twiki/bin/view/Atlas/AFPRunning#Tracker
  if (hitLink <=3 ) {	// not checking if larger or equal 0, because this is uint, which is always positive
    // set near station ID selecting side based on ROB
    if (robID == AFP_ROBID::sideA)
      xAODSiHit->setStationID(xAOD::AFPStationID::nearA);
    else if (robID == AFP_ROBID::sideC || robID == AFP_ROBID::sideC_2016)
      xAODSiHit->setStationID(xAOD::AFPStationID::nearC);
    else  {
      ATH_MSG_WARNING("Unrecognised robID: in dec="<<std::dec<<robID<<",  in hex=0x"<<std::hex<<robID<<std::dec);
      xAODSiHit->setStationID(-1);
    }
  } else if (hitLink >= 8 && hitLink <= 11) {
    // set far station ID selecting side based on ROB
    if (robID == AFP_ROBID::sideA)
      xAODSiHit->setStationID(xAOD::AFPStationID::farA);
    else if (robID == AFP_ROBID::sideC || robID == AFP_ROBID::sideC_2016)
      xAODSiHit->setStationID(xAOD::AFPStationID::farC);
    else  {
      ATH_MSG_WARNING("Unrecognised robID: in dec="<<std::dec<<robID<<",  in hex=0x"<<std::hex<<robID<<std::dec);
      xAODSiHit->setStationID(-1);
    }
  } else {
    ATH_MSG_WARNING("Unrecognised value in hitLink = " << hitLink);
    xAODSiHit->setStationID(-1);
  } // end of switch
  
  xAODSiHit->setPixelLayerID( hitLink%4 );
  xAODSiHit->setPixelColIDChip( data.column() );
  xAODSiHit->setPixelRowIDChip( data.row() );

  const unsigned int ToT_value = decodeTimeOverThresholdSi( data.timeOverThreshold(), data.hitDiscConfig() );
  xAODSiHit->setTimeOverThreshold(ToT_value);
  xAODSiHit->setDepositedCharge( m_totToChargeTransformation(ToT_value) );

  ATH_MSG_DEBUG("have AFPSiHit: pixelLayerID "<<xAODSiHit->pixelLayerID()<<", pixelColIDChip "<<xAODSiHit->pixelColIDChip()<<", pixelRowIDChip "<<xAODSiHit->pixelRowIDChip()<<", timeOverThreshold "<<xAODSiHit->timeOverThreshold()<<", depositedCharge "<<xAODSiHit->depositedCharge());
}

void AFP_Raw2DigiTool::setBarAndTrainID(xAOD::AFPToFHit* tofHit, const EventContext& ctx) const
{
  // mapping is implemented according to https://twiki.cern.ch/twiki/bin/view/Atlas/AFPHPTDC#Channel_Mapping

  const int hptdcChannel = tofHit->hptdcChannel();
  const int hptdcID = tofHit->hptdcID();
  
  if((hptdcChannel<0 || hptdcChannel>=12) || (hptdcID<=0 || hptdcID>=3))
  {
  	ATH_MSG_WARNING("Invalid hptdcID "<<hptdcID<<", or hptdcChannel "<<hptdcChannel);
  	return;
  }
  
  int RunID = (ctx.eventID().run_number() > 370000 ? (ctx.eventID().run_number() > 436657 ? 2 : 1) : 0);
  
  int trainID=m_channel2train[RunID][hptdcID-1][hptdcChannel];
  int barID=m_channel2bar[RunID][hptdcID-1][hptdcChannel];
  
  if(trainID>=0 && barID>=0)
  {
    tofHit->setTrainID(trainID);
    tofHit->setBarInTrainID(barID);
  }
  else if(trainID==-2 && barID==-2)
  {
    // this is known problem for Run2, channels 1 and 4
    ATH_MSG_DEBUG("Unrecognised hptdcChannel "<<hptdcChannel<<", hptdcID "<<hptdcID);
  }
  else if(trainID==-1 && barID==-1)
  {
    // this is not known
    ATH_MSG_WARNING("Unrecognised hptdcChannel "<<hptdcChannel<<", hptdcID "<<hptdcID);
  }
  else
  {
    // this shouldn't happen
    ATH_MSG_WARNING("Unrecognised trainID "<<trainID<<" or barID "<<barID<<"; run "<<ctx.eventID().run_number()<<", hptdcChannel "<<hptdcChannel<<", hptdcID "<<hptdcID);
  }
}

unsigned int AFP_Raw2DigiTool::decodeTimeOverThresholdSi (const unsigned int input, const unsigned int discConfig) const
{
  if (input <= 13)
    return input + discConfig + 1;
  else if (input == 14) 
    return discConfig;
  else
    return 0;
}
