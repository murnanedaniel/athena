/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration 
*/
#include "InDetSecVxFinderTool/JetFitterTrackSelectorTool.h"
#include <cassert>

using namespace InDet;

    JetFitterTrackSelectorTool::JetFitterTrackSelectorTool(const std::string &t, const std::string &n, const IInterface *p)
            : AthAlgTool(t, n, p)
    {
        declareInterface< JetFitterTrackSelectorTool >(this);
    }

    JetFitterTrackSelectorTool::~JetFitterTrackSelectorTool() = default;

    StatusCode JetFitterTrackSelectorTool::initialize() {

	if ( m_trkFilter.retrieve().isFailure() ) {
	  msg(MSG::ERROR) << " Unable to retrieve InDet::InDetDetailedTrackSelectorTool" << endmsg;
	  return StatusCode::FAILURE;
	}

	if ( m_jetFitterUtils.retrieve().isFailure() ) {
	  msg(MSG::ERROR) << " Unable to retrieve InDet::InDetJetFitterUtils/InDetJetFitterUtils" << endmsg;
	  return StatusCode::FAILURE;
	}

	if ( m_extrapolator.retrieve().isFailure() ) {
	  msg(MSG::ERROR) << " Unable to retrieve Trk::Extrapolator/InDetExtrapolator" << endmsg;
          return StatusCode::FAILURE;
	}

        return StatusCode::SUCCESS;
    }


    const Trk::SelectedTracksInJet* JetFitterTrackSelectorTool::doTrackSelection( const xAOD::Vertex &primaryVertex,
										  const TLorentzVector &jetMomentum,
										  const std::vector<const xAOD::IParticle *> &inputTracks) const {
      // perform the track selection
      // step 1, apply a track filter like "InDet::InDetDetailedTrackSelectorTool"
      // step 2, calculate the compatibility of filtered tracks with primary vertex
      // use this to deduce primaryTracks and secondaryTracks

      ATH_MSG_DEBUG( "Doing track selection on " << inputTracks.size() << " tracks ... " );

      // We need to use normal pointers instead of smart pointers since the code breaks.
      // We have to fix this issue in the future
      //      if ( m_selectedTracks != nullptr ) delete m_selectedTracks; // May this break the code?
      Trk::SelectedTracksInJet *selectedTracks = new Trk::SelectedTracksInJet();

      // Vectors of Trk::ITrackLink to be given to m_selectedTracks once we understand if they are primary of secondary tracks
      std::vector< const Trk::ITrackLink* > primaryTrackLinks;
      std::vector< const Trk::ITrackLink* > secondaryTrackLinks;

      // Running on input tracks
      std::vector<const xAOD::IParticle *>::const_iterator trk_iter = inputTracks.begin();
      std::vector<const xAOD::IParticle*>::const_iterator   trk_end = inputTracks.end();  

      int counter = 0;
      for ( ; trk_iter != trk_end; ++trk_iter ) {
	// Convert xAOD::IParticle to xAOD::TrackParticle
	const xAOD::TrackParticle * tmp = dynamic_cast< const xAOD::TrackParticle* > ( *trk_iter );
	assert( tmp != nullptr ); // in principle should really check that inputTracks only contains TrackParticle objects

	// Compute compatibility and understand track type
	// -1: track filter failed
	//  0: extrapolation of MeasuredPerigee failed
	//  1: primary
	//  2: secondary
	int type = computeTrackCompatibility( primaryVertex,jetMomentum,*tmp );

	// Create Trk::ITrackLink collections to be given to selected tracks
	if (type==1 || type==2) {
           ElementLink< xAOD::TrackParticleContainer > linkTP;
           linkTP.setElement( tmp );
           Trk::LinkToXAODTrackParticle* link= new Trk::LinkToXAODTrackParticle( linkTP );

           if ( type == 1) primaryTrackLinks.push_back( link );
           else if ( type == 2 ) secondaryTrackLinks.push_back( link );
	}
	else {
           continue;
	}

	// How many tracks we are selecting
	counter++;
      }

      ATH_MSG_DEBUG( "  Total of selected tracks: "<< counter );

      selectedTracks->setPrimaryTrackLinks( primaryTrackLinks );
      selectedTracks->setSecondaryTrackLinks( secondaryTrackLinks );
      return selectedTracks;
   }

    int JetFitterTrackSelectorTool::computeTrackCompatibility( const xAOD::Vertex &primaryVertex,
							    const TLorentzVector &jetMomentum,
							       const xAOD::TrackParticle &track ) const {

      // Decorators for tracks
      std::string toolname = this->name();
      std::string delimiter = "_";
      std::string::size_type firstDelimiter = toolname.find(delimiter);
      std::string sub = toolname.substr(0, firstDelimiter);
      std::string decoratorName = std::string("JetFitter_TrackCompatibility_") + sub;
      SG::AuxElement::Decorator< float > compatibilityDecorator(decoratorName.c_str());
      
      // Apply track filter
      if ( !m_trkFilter->decision( track, &primaryVertex ) ) {
        compatibilityDecorator ( track ) = 0.;
        return -1;
      }

      // Recomputing Perigee w.r.t PV
      Trk::PerigeeSurface mySurface( primaryVertex.position() );
      std::unique_ptr<const Trk::TrackParameters>  myMeasuredPerigee(m_extrapolator->extrapolate(
          Gaudi::Hive::currentContext(),track,mySurface ));

      if ( !myMeasuredPerigee) {
        ATH_MSG_DEBUG( " Extrapolation to primary vertex failed. Skipping track " );
        compatibilityDecorator ( track ) = 0.;
        return 0;
      }
      
      
      // Prepare for using jetFitterUtils (for the computation of the compatibility)
      // Is this conrvertion really necessary?
      Trk::RecVertex primaryVertexRecVertex( primaryVertex.position(), 
					     primaryVertex.covariancePosition(), 
					     primaryVertex.numberDoF(), 
					     primaryVertex.chiSquared());
      
      Amg::Vector3D jetMomSpatial( jetMomentum.X(),jetMomentum.Y(),jetMomentum.Z() );
      double compatibilityValue = m_jetFitterUtils->compatibility( *myMeasuredPerigee,primaryVertexRecVertex ).first;
      compatibilityValue = fabs( compatibilityValue ) * m_jetFitterUtils->get3DLifetimeSignOfTrack( *myMeasuredPerigee,
												    jetMomSpatial,
												    primaryVertexRecVertex );
      
      // Decorate
      ATH_MSG_DEBUG( "compatibilityValue = " << compatibilityValue );
      compatibilityDecorator ( track ) = compatibilityValue;

      // Understand if primary or secondary track particle
      double cutCompatibilityPVforPosTracks = m_cutCompatibilityPrimaryVertexForPositiveLifetimeTracks;
      double cutCompatibilityPVforNegTracks = m_cutCompatibilityPrimaryVertexForNegativeLifetimeTracks;
      
      if ( m_revertFromPositiveToNegativeTags ) {
	cutCompatibilityPVforNegTracks = m_cutCompatibilityPrimaryVertexForPositiveLifetimeTracks;
	cutCompatibilityPVforPosTracks = m_cutCompatibilityPrimaryVertexForNegativeLifetimeTracks;
      }
      
      if ( ( compatibilityValue < 0 && 
	     TMath::Prob( fabs( compatibilityValue ),2 ) < cutCompatibilityPVforNegTracks) ||
	   ( compatibilityValue >= 0 && 
	     TMath::Prob( fabs( compatibilityValue ),2 ) < cutCompatibilityPVforPosTracks ) ) 
	return 2;
      else return 1;

    }



