/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/// @file   AFPTDBasicTool.cxx
/// @author Leszek Adamczyk <leszek.adamczyk@agh.edu.pl>
/// @date   2021-11-11
/// 
/// @brief  Implementation file for AFPTDBasicTool used in tracks reconstruction.


// FrameWork includes
#include "GaudiKernel/IToolSvc.h"
#include "CLHEP/Units/SystemOfUnits.h"

// AFP_LocReco includes
#include "AFP_LocReco/AFPTDBasicTool.h"

AFPTDBasicTool::AFPTDBasicTool( const std::string& type, 
                                const std::string& name, 
                                const IInterface* parent ) : 
  base_class ( type, name, parent )
{

}

StatusCode AFPTDBasicTool::initialize()
{
  // print information about initialised stations
  ATH_MSG_DEBUG("AFPTDBasicTool::initialize()");
  ATH_MSG_INFO("Station with ID="<<m_stationID <<" will have minimum number of "<<m_minHitsNumber <<" bars.");
  ATH_MSG_INFO("Maximal length of signal  at which bar can be joined to the track  m_maxAllowedLength = "<<m_maxAllowedLength);
  
  CHECK( m_tofLocParamDBTool.retrieve() );
  CHECK( m_hitContainerKey.initialize( SG::AllowEmpty ) );

  return StatusCode::SUCCESS;
}


StatusCode AFPTDBasicTool::finalize()
{
  ATH_MSG_DEBUG("AFPTDBasicTool::finalize()");
  return StatusCode::SUCCESS;
}


void AFPTDBasicTool::fillTrainWithBars(std::vector<const xAOD::AFPToFHit*> my_trainBars[4], SG::ReadHandle<xAOD::AFPToFHitContainer>& hitContainer) const
{
  ATH_MSG_DEBUG("AFPTDBasicTool::fillTrainWithBars");
  ATH_MSG_DEBUG("Total number of AFP ToF hits read in = " << hitContainer->size());
  // retrieve ToF bars
  try {
    // fill station with ToF hits
    for (const xAOD::AFPToFHit* theHit : *hitContainer)
      if (theHit->stationID() == m_stationID) // check if hit is from the correct station
        my_trainBars[theHit->trainID()].push_back(theHit);
  }
  catch (const std::out_of_range& outOfRange) {
    ATH_MSG_WARNING("Bar with station outside expected range. Aborting track reconstruction.");
    clearTrains(my_trainBars);
  }
}

StatusCode AFPTDBasicTool::reconstructTracks(std::unique_ptr<xAOD::AFPToFTrackContainer>& outputContainer, const EventContext& ctx) const
{

  if (m_hitContainerKey.empty()) {
    // this is allowed, there might be no AFP data in the input
    return StatusCode::SUCCESS;
  }

  SG::ReadHandle<xAOD::AFPToFHitContainer> hitContainer( m_hitContainerKey, ctx );
  if(!hitContainer.isValid())
  {
      // this is allowed, there might be no AFP data in the input
      return StatusCode::SUCCESS;
  }
  
  // prepare list for storing  reconstructed ToF tracks
  std::list<AFPTDBasicToolTrack> reconstructedTracks;

  std::vector<const xAOD::AFPToFHit*> my_trainBars[4];

  clearTrains(my_trainBars);
  fillTrainWithBars(my_trainBars, hitContainer);

  ATH_MSG_DEBUG( "Number of AFP ToF hits in each train = " << my_trainBars[0].size()<<", "<<my_trainBars[1].size()<<", "<<my_trainBars[2].size()<<", "<<my_trainBars[3].size());

  // ===== do tracks reconstruction =====
  
  nlohmann::json dataTLP;
  bool dataTLP_init{false};
  
  for (int k=0; k<4; k++)
  {
    unsigned int TrSize = my_trainBars[k].size();
    double TrTime = 0.;
    double weight = 0;
    unsigned int TrSat = 0;
    if ( TrSize >= m_minHitsNumber )
    {
      for (unsigned int l=0; l<TrSize; l++)
      {
        // for l-th hit in the k-th train, not necessary l-th bar
        if ( my_trainBars[k].at(l)->pulseLength() < m_maxAllowedLength) 
        {
          if(!dataTLP_init)
          {
            // read from DB only if necessary
            dataTLP=m_tofLocParamDBTool->parametersData(ctx);
            dataTLP_init=true;
          }
          
          const AFP::ToFLocParamData TLP=m_tofLocParamDBTool->parameters(dataTLP, m_stationID, my_trainBars[k].at(l)->trainID(), my_trainBars[k].at(l)->barInTrainID());
          
          TrTime += (my_trainBars[k].at(l)->time()-TLP.barTimeOffset())*TLP.barWeight();
          weight += TLP.barWeight();
        }
        else
        {
          TrSat++;
        }
      }
    // time average
    if( TrSize!=TrSat) TrTime /= weight;
    
    ATH_MSG_DEBUG("Track reconstruction complete: stationID = " + std::to_string(m_stationID) + ", train time = " + std::to_string(TrTime) + ", train size = " + std::to_string(TrSize));
    reconstructedTracks.emplace_back(m_stationID,k,TrTime, TrSize, TrSat);
    AFPTDBasicToolTrack& theTrack = reconstructedTracks.back();
    for(unsigned int l=0; l<TrSize; l++) theTrack.addBar(my_trainBars[k].at(l));
    }
  }

  // === Save result to xAOD ===
  // save to xAOD

  for (const AFPTDBasicToolTrack& track : reconstructedTracks)
    saveToXAOD(track, outputContainer, hitContainer);  
  return StatusCode::SUCCESS;
}

void AFPTDBasicTool::saveToXAOD (const AFPTDBasicToolTrack& recoTrack, std::unique_ptr<xAOD::AFPToFTrackContainer>& containerToFill, SG::ReadHandle<xAOD::AFPToFHitContainer>& hitContainer) const
{
    ATH_MSG_DEBUG("AFPTDBasicTool::saveToXAOD");
    auto track = containerToFill->push_back(std::make_unique<xAOD::AFPToFTrack>());
  
    track->setStationID(recoTrack.m_stationID);
    track->setTrainID(recoTrack.m_trainID);
    track->setTrainTime(recoTrack.m_trainTime);
    track->setTrainSize(recoTrack.m_trainSize);
    track->setTrainNSat(recoTrack.m_trainNSat);
    track->setAlgID(0);
// add links to bars
    ATH_MSG_DEBUG("Track time: (time="<<track->trainTime()<<", size="<<track->trainSize()<<", train id="<<track->trainID()<<")   station: (st="<<track->stationID()<<")");
    for (const xAOD::AFPToFHit* theHit : recoTrack.barInTrain()) {
        ElementLink< xAOD::AFPToFHitContainer >* barLink = new ElementLink< xAOD::AFPToFHitContainer >; // will be taken over by the xAODTrain and deleted
        barLink->toContainedElement(*hitContainer, theHit);
        track->addBar(*barLink);

        ATH_MSG_DEBUG("bar time="<<theHit->time()<<" bar in train ="<<theHit->barInTrainID()<<", trin id ="<<theHit->trainID()<<")");
    }

}

