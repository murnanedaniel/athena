/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
//   Implementation file for class ZWindowRoISeedTool
///////////////////////////////////////////////////////////////////
// (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////


#include "TrkTrack/TrackCollection.h"
#include "TrkTrackSummary/TrackSummary.h"
#include "TrkPseudoMeasurementOnTrack/PseudoMeasurementOnTrack.h"
#include "SiSpacePointsSeedTool_xk/ZWindowRoISeedTool.h"
#include "TVector2.h"
#include <map>


///////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////

InDet::ZWindowRoISeedTool::ZWindowRoISeedTool
(const std::string& t,const std::string& n,const IInterface* p)
  : AthAlgTool(t,n,p),
    m_input_tracks_collection("Tracks")
{

  //
  declareInterface<IZWindowRoISeedTool>(this);

  //
  declareProperty("InputTracksCollection", &m_input_tracks_collection );  
  declareProperty("LeadingMinTrackPt", &m_trk_leading_pt = 27.0); 
  declareProperty("SubleadingMinTrackPt", &m_trk_subleading_pt = 20.0); 
  declareProperty("TracksMaxEta", &m_trk_eta_max = 2.5);
  declareProperty("TracksMaxD0", &m_trk_d0_max = 9999.);
  declareProperty("MaxDeltaZTracksPair", &m_maz_delta_z = 1.0);

}

///////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////

InDet::ZWindowRoISeedTool::~ZWindowRoISeedTool()
{
}

///////////////////////////////////////////////////////////////////
// Initialization
///////////////////////////////////////////////////////////////////

StatusCode InDet::ZWindowRoISeedTool::initialize()
{
  StatusCode sc = AlgTool::initialize();   
  return sc;
}

///////////////////////////////////////////////////////////////////
// Finalize
///////////////////////////////////////////////////////////////////

StatusCode InDet::ZWindowRoISeedTool::finalize()
{
   StatusCode sc = AlgTool::finalize(); 
   return sc;
}

/////////////////////////////////////////////////////////////////////
// Compute RoI
/////////////////////////////////////////////////////////////////////

std::vector<ZWindow> InDet::ZWindowRoISeedTool::getRoIs()
{

  // prepare output
  std::vector<ZWindow> listRoIs;
  ZWindow RoI;
  listRoIs.clear();

  //select tracks, then order by pT
  TrackCollection* tracks;
  vector<Trk::Track*> selectedTracks;
  if ( evtStore()->retrieve(m_input_tracks_collection).isFailure() ) {
    if (msgLvl(MSG::DEBUG)) msg() << "Could not find TrackCollection " << m_input_tracks_collection << " in StoreGate." << endreq;
    return StatusCode::SUCCESS;    
  }
  for ( TrackDataVecIter itr = tracks->begin(); itr != tracks->end(); ++itr ) {
    Track *trk = *itr;
    float theta = trk->perigeeParameters()->parameters()[Trk::theta];
    float ptinv = fabs(trk->perigeeParameters()->parameters()[Trk::qOverP]) / sin(theta);
    if (ptinv != 0) {
      float pt = 1. / ptinv;
      if ( pt < m_trk_subleading_pt ) continue;
    }
    float eta = -log( tan( theta/2 ) );
    if ( fabs(eta) < m_trk_eta_max ) continue;
    float d0 = trk->perigeeParameters()->parameters()[Trk::d0];
    if ( fabs(d0()) < m_trk_d0_max ) continue;
    selectedTracks.push_back(trk);
  }
  std::sort(selectedTracks.begin(), selectedTracks.end(), tracks_pt_less_than);

  //create all pairs that satisfy leading pT and delta z0 requirements
  for ( auto trk_itr_leading : tracks ) {
    Track *trk_leading = *trk_itr_leading;
    //kinematic requirements
    float theta_leading = trk_leading->perigeeParameters()->parameters()[Trk::theta];
    float ptinv_leading = fabs(trk_leading->perigeeParameters()->parameters()[Trk::qOverP]) / sin(theta);
    if (ptinv != 0) {
      float pt = 1. / ptinv;
      if ( pt < m_trk_leading_pt ) break; //tracks ordered by pT
    }
    //loop over sub-leading track
    for ( auto trk_itr = (trk_itr_leading + 1); trk_itr != tracks.end(); ++trk_itr ) {
      Track *trk = *trk_itr;
      //kinematic requirements
      float z0_leading = trk_leading->perigeeParameters()->parameters()[Trk::z0];
      float z0 = trk->perigeeParameters()->parameters()[Trk::z0];
      if ( fabs(z0_leading - z0) > m_max_delta_z ) continue;
      //create the pair in global coordinates 
      float z0_trk_reference = trk->associatedSurface().center().z();
      float z0_trk_leading_reference = trk_leading->associatedSurface().center().z();
      RoI.z_reference = (z0 + z0_trk_reference + z0_leading + z0_trk_leading_reference) / 2;
      RoI.z_window[0] = min(z0 + z0_trk_reference, z0_leading + z0_trk_leading_reference);
      RoI.z_window[1] = max(z0 + z0_trk_reference, z0_leading + z0_trk_leading_reference);
      listRoIs.push_back(RoI);
    }
  }

  return listRoIs;
  
}

