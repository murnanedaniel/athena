/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file InDetPhysHitDecoratorTool.cxx
 * @author shaun roe
**/

#include "InDetPhysHitDecoratorTool.h"
#include "xAODTracking/TrackParticle.h"
//#include "GeneratorUtils/PIDUtils.h"
#include "TrkParameters/TrackParameters.h" //Contains typedef to Trk::CurvilinearParameters
#include "TrkToolInterfaces/ITrackHoleSearchTool.h"
#include "TrkToolInterfaces/IUpdator.h"
#include "InDetRIO_OnTrack/SiClusterOnTrack.h"
#include "TrkEventPrimitives/ResidualPull.h"
#include "TrkToolInterfaces/IResidualPullCalculator.h"
#include "TrkTrack/TrackCollection.h"
//for the identifiers
#include "AtlasDetDescr/AtlasDetectorID.h"
#include "InDetIdentifier/PixelID.h"
#include "InDetIdentifier/SCT_ID.h"
#include "InDetIdentifier/TRT_ID.h"
#include "InDetPrepRawData/SiCluster.h"
//
#include <tuple>
#include <limits>

//ref: ​https://svnweb.cern.ch/trac/atlasoff/browser/Tracking/TrkEvent/TrkParametersBase/trunk/TrkParametersBase/CurvilinearParametersT.h


InDetPhysHitDecoratorTool::InDetPhysHitDecoratorTool(const std::string& type,const std::string& name,const IInterface* parent):
AthAlgTool(type,name,parent), 
m_holeSearchTool("InDet::InDetTrackHoleSearchTool"),
m_updatorHandle("Trk::KalmanUpdator/TrkKalmanUpdator"),
m_residualPullCalculator("Trk::ResidualPullCalculator/ResidualPullCalculator"),
m_ptThreshold(0.8), m_isUnbiased(false), m_doUpgrade(false){
declareInterface<IInDetPhysValDecoratorTool>(this);

declareProperty("InDetTrackHoleSearchTool"     , m_holeSearchTool);
declareProperty("Updator"     , m_updatorHandle);
declareProperty("ResidualPullCalculator", m_residualPullCalculator);
  //do I need to retrieve the 'Tracks' container?
}

InDetPhysHitDecoratorTool::~InDetPhysHitDecoratorTool (){
//nop
}

StatusCode 
InDetPhysHitDecoratorTool::initialize(){
	StatusCode sc = AlgTool::initialize();
  if (sc.isFailure()) return sc;
  ATH_CHECK( m_holeSearchTool.retrieve() );
  if ( not (m_updatorHandle.empty())) {ATH_CHECK(m_updatorHandle.retrieve());}
  if ( not (m_holeSearchTool.empty())) {ATH_CHECK(m_holeSearchTool.retrieve());}
  
   //ID Helper
  m_idHelper = new AtlasDetectorID;

  // Get the dictionary manager from the detector store
  ATH_CHECK(detStore()->retrieve(m_idHelper, "AtlasID"));
  ATH_CHECK(detStore()->retrieve(m_pixelID, "PixelID"));
  ATH_CHECK(detStore()->retrieve(m_sctID, "SCT_ID"));
  
  if (!m_doUpgrade) {ATH_CHECK(detStore()->retrieve(m_trtID, "TRT_ID"));}
  if (m_residualPullCalculator.empty()) {
    ATH_MSG_INFO("No residual/pull calculator for general hit residuals configured.");
    ATH_MSG_INFO("It is recommended to give R/P calculators to the det-specific tool handle lists then.");
  } else if (m_residualPullCalculator.retrieve().isFailure()) {
    msg(MSG::FATAL) << "Could not retrieve "<< m_residualPullCalculator<<" (to calculate residuals and pulls) "<< endreq;
    return StatusCode::FAILURE;
  } else {
    ATH_MSG_INFO("Generic hit residuals & pulls will be calculated in one or both available local coordinates");
  }
	return sc;
}

StatusCode 
InDetPhysHitDecoratorTool::finalize  (){
	return StatusCode::SUCCESS;
}

bool
InDetPhysHitDecoratorTool::decorateTrack(const xAOD::TrackParticle & particle, const std::string& prefix){
	typedef std::tuple<int,int,int, float, float, float, float, int> SingleResult_t;
	typedef std::vector<SingleResult_t> TrackResult_t;
	const float invalidFloat(std::numeric_limits<float>::quiet_NaN());
	//const float invalidDouble(std::numeric_limits<double>::quiet_NaN());
	const float invalidRes(invalidFloat), invalidPull(invalidFloat);
	const int invalidDetector(-1);
	const int invalidRegion(-1);
	const int invalidLayer(-1);
	const int invalidWidth(-1);
	const SingleResult_t invalidResult=std::make_tuple(invalidDetector,invalidRegion,invalidLayer, invalidRes, invalidPull,invalidRes,invalidPull,invalidWidth);
  //get element link to the original track
  const ElementLink< TrackCollection >& trackLink = particle.trackLink();//using xAODTracking-00-03-09, interface has changed later
  if (trackLink.isValid()){
    ATH_MSG_VERBOSE ("Track link found " );
    const double pt = particle.pt(); 
    if (pt > m_ptThreshold){
    ATH_MSG_VERBOSE ("pt is over threshold " );
			std::unique_ptr<const Trk::Track> trackWithHoles(m_holeSearchTool->getTrackWithHoles(**trackLink));
			const int numberOfHits((trackWithHoles->trackStateOnSurfaces())->size());
			//ATH_MSG_INFO ("number of Hits "<<numberOfHits );
			TrackResult_t result; result.reserve(numberOfHits);
			for (auto &thisTrackState: *(trackWithHoles->trackStateOnSurfaces())){
				SingleResult_t thisResult(invalidResult);
				if (not thisTrackState) continue;
				const Trk::MeasurementBase* mesb=thisTrackState->measurementOnTrack();
				if (not mesb) {
				  ATH_MSG_DEBUG ("intermediate mesb is NULL");
				} else {
				  ATH_MSG_VERBOSE ("intermediate mesb is ok");
				}
				const Trk::RIO_OnTrack* hit = mesb ? dynamic_cast<const Trk::RIO_OnTrack*>(mesb) : 0;
				const Trk::TrackParameters* biasedTrackParameters = thisTrackState->trackParameters();
				if (biasedTrackParameters) ATH_MSG_VERBOSE ("biased track parameters ok");
				Identifier surfaceID;
				if (mesb && mesb->associatedSurface().associatedDetectorElement()) {
          surfaceID = mesb->associatedSurface().associatedDetectorElement()->identify();
          ATH_MSG_DEBUG("Surface ID found");
        } else { // holes, perigee
          if (not biasedTrackParameters ) {
            msg(MSG::DEBUG) << "The track parameters are not valid for this track state (the pointer is null)" << endreq;
            continue;
          }
          surfaceID = biasedTrackParameters->associatedSurface().associatedDetectorElementIdentifier();
          ATH_MSG_VERBOSE("Surface ID found for holes etc");
        }
				ATH_MSG_VERBOSE("checking mesb and track parameters");
				if (mesb && biasedTrackParameters) {
					ATH_MSG_DEBUG("mesb and biased track parameters are ok");
					const Trk::TrackParameters *trackParameters = (! thisTrackState->type(Trk::TrackStateOnSurface::Outlier)) ? getUnbiasedTrackParameters(biasedTrackParameters,mesb) : biasedTrackParameters;
				  if (not trackParameters){
				    ATH_MSG_DEBUG("unbiased track parameters pointer is NULL");
				  }
				  Trk::ResidualPull::ResidualType resType = (m_isUnbiased) ? (Trk::ResidualPull::Unbiased):(Trk::ResidualPull::Biased);
				  std::unique_ptr<const Trk::ResidualPull> residualPull(m_residualPullCalculator->residualPull(hit,trackParameters,resType));
				  ATH_MSG_VERBOSE("checking residual pull");
				  if (not residualPull){
				    ATH_MSG_DEBUG("residualPull is NULL");
				    continue;
				  }
				  float residualLocX = 1000.*residualPull->residual()[Trk::loc1]; // residuals in microns
				  float pullLocX = residualPull->pull()[Trk::loc1];
				  float residualLocY(invalidFloat),  pullLocY(invalidFloat);//NaN by default
				  if (residualPull->dimension() > 1){
				  	residualLocY = 1000.*residualPull->residual()[Trk::loc2];
				  	pullLocY = residualPull->pull()[Trk::loc2];
				  }
				  Subdetector det(INVALID_DETECTOR);
				  Region r(INVALID_REGION);
				  int iLayer(invalidLayer);
				  const bool successfulIdentification = decideDetectorRegion(surfaceID, det, r, iLayer);
				  if (not successfulIdentification) {
				    ATH_MSG_DEBUG("Could not identify surface");
				  	continue;
				  }
				  //int width = 1; //check original code
				  int phiWidth(-1);
				  //int zWidth(-1);
					//copy-paste from original
					if (hit && m_isUnbiased) {
						// Cluster width determination
						if((det == PIXEL) or  (det==SCT)) {
							const InDet::SiCluster* pCluster = dynamic_cast <const InDet::SiCluster*>(hit->prepRawData());
							if(pCluster){
								InDet::SiWidth width = pCluster->width();
								phiWidth = int(width.colRow().x());
								//zWidth = int(width.colRow().y());
							}
						}
					}
					//end copy-paste
				  thisResult=std::make_tuple(det, r, iLayer, residualLocX, pullLocX, residualLocY, pullLocY, phiWidth);
				  ATH_MSG_INFO ("**result "<<iLayer<<", "<<residualLocX<<", "<<pullLocX<<", "<<residualLocY<<", "<<pullLocY<<", "<<phiWidth );
				  result.push_back(thisResult);
				  //must delete the pointers?
				} else {
					if (not mesb) ATH_MSG_INFO("mesb not ok");
					if (not biasedTrackParameters) ATH_MSG_INFO("biasedTrackParameters were not found");
				}
			}//end of for loop
			if (not result.empty()){
			  particle.auxdecor<TrackResult_t>(prefix+"hitResiduals") = result;
			  return true;
			} else {
			  ATH_MSG_INFO("No hit residual added");
			}
		} 
	} else {
	  ATH_MSG_DEBUG ("No valid track link found " );
	} 
	return false;
}

bool
InDetPhysHitDecoratorTool::decideDetectorRegion(const Identifier & id, Subdetector & det, Region & r, int & layer){
	bool success(false);
	const int normalBarrel(0);
	const int upgradedBarrel(1);
	const int normalTrtBarrel(1);
	const int dbm(2);

	det = INVALID_DETECTOR;//default
	r = INVALID_REGION;
	int bec(-100);
	
	//following the logic in the original code, should be reviewed!
	if (m_idHelper->is_pixel(id)) {
	  bec = abs(m_pixelID->barrel_ec(id));
	  if (bec == dbm)
	    det=DBM;
	  else
	    det=PIXEL;
	}
	if (m_idHelper->is_sct(id)) det=SCT;
	if (not m_doUpgrade and m_idHelper->is_trt(id)) det=TRT;
	//
	//check this specifically
	if (det==PIXEL) {
		bec = abs(m_pixelID->barrel_ec(id));
		r= (bec==normalBarrel)?(BARREL):(ENDCAP);
		layer = m_pixelID->layer_disk(id);
		if (layer == 0) det=BLAYER;
	}
	/** cf. Miriam's code 
	if (det==PIXEL) {
	  r= (bec==normalBarrel)?(BARREL):(ENDCAP);
	  if (m_pixelID->layer_disk(id) == 0) det=BLAYER;
	}
	**/
	if (det==DBM){
	  r= (bec < 0) ? (BARREL) : (ENDCAP) ;
	}
	if (det==SCT) {
		bec = abs(m_sctID->barrel_ec(id));
		if (not m_doUpgrade){
		  r = (bec==normalBarrel)?(BARREL):(ENDCAP);
		} else {
			r = (bec==upgradedBarrel)?(BARREL):(ENDCAP);
		}
		layer = m_sctID->layer_disk(id);
	}
	if (det==TRT) {
		bec = abs(m_trtID->barrel_ec(id));
		r = (bec==normalTrtBarrel)?(BARREL):(ENDCAP);
		layer = m_trtID->layer_or_wheel(id);
	}
	success = (det != INVALID_DETECTOR) and (r != INVALID_REGION);
	
	return success;
}

const Trk::TrackParameters*
InDetPhysHitDecoratorTool::getUnbiasedTrackParameters(const Trk::TrackParameters* trkParameters, const Trk::MeasurementBase* measurement ) {
	static bool alreadyWarned(false);
	const Trk::TrackParameters* unbiasedTrkParameters(0);
	if (!m_updatorHandle.empty() && (m_isUnbiased) ) {
    if ( trkParameters->covariance() ) {
      // Get unbiased state
      unbiasedTrkParameters =m_updatorHandle->removeFromState( *trkParameters,measurement->localParameters(),measurement->localCovariance());
      if (!unbiasedTrkParameters) {
    		msg(MSG::WARNING) << "Could not get unbiased track parameters, use normal parameters" << endreq;
    		m_isUnbiased = false;
      }
    } else if (not alreadyWarned) {
      // warn only once!
      msg(MSG::WARNING) << "TrackParameters contain no covariance, unbiased track states can not be calculated (ie. pulls and residuals will be too small)" << endreq;
      alreadyWarned = true;
      m_isUnbiased = false;
    } else {
      m_isUnbiased = false;
    }// end if no measured track parameter
  }
  return unbiasedTrkParameters;                                 
}

  /**
for (; TSOSItr != trackWithHoles->trackStateOnSurfaces()->end(); ++TSOSItr) {
	const Trk::MeasurementBase* mesb=(*TSOSItr)->measurementOnTrack();
        const Trk::RIO_OnTrack* hit = mesb ? dynamic_cast<const Trk::RIO_OnTrack*>(mesb) : 0;
	const Trk::TrackParameters* biasedTrackParameters = (*TSOSItr)->trackParameters();
                
	if (mesb && biasedTrackParameters) {
		const Trk::TrackParameters *trackParameters = (!(*TSOSItr)->type(Trk::TrackStateOnSurface::Outlier)) ?getUnbiasedTrackParameters(biasedTrackParameters,mesb) : biasedTrackParameters;
	
		Trk::ResidualPull::ResidualType resType = (m_isUnbiased) ? (Trk::ResidualPull::Unbiased):(Trk::ResidualPull::Biased);
		const auto_ptr<const Trk::ResidualPull> residualPull(m_residualPullCalculator->residualPull(hit,trackParameters,resType));
		residualLocX = 1000*residualPull->residual()[Trk::loc1]; // residuals in microns
		m_residualx_pixel_barrel->Fill(residualLocX);
	}
}
}      
        **/
    
  





/**
const Trk::TrackParameters*
IDStandardPerformance::getUnbiasedTrackParameters(const Trk::TrackParameters* trkParameters
                                                , const Trk::MeasurementBase* measurement ) {
  const Trk::TrackParameters *unbiasedTrkParameters = 0;

  // ------------------------------------
  // try to get measured track parameters
  // ------------------------------------

  //const Trk::MeasuredTrackParameters *measuredTrkParameters =
  //dynamic_cast<const Trk::MeasuredTrackParameters*>(trkParameters);

  if (!m_updatorHandle.empty() && (m_isUnbiased==1) ) {
    if ( trkParameters->covariance() ) {
      // Get unbiased state
      unbiasedTrkParameters =
    m_updatorHandle->removeFromState( *trkParameters,
                    measurement->localParameters(),
                    measurement->localCovariance());

      if (!unbiasedTrkParameters) {
    msg(MSG::WARNING) << "Could not get unbiased track parameters, "
        <<"use normal parameters" << endreq;
    m_isUnbiased = 0;
      }
    } else if (!m_UpdatorWarning) {
      // warn only once!
      msg(MSG::WARNING) << "TrackParameters contain no covariance: "
      <<"Unbiased track states can not be calculated "
      <<"(ie. pulls and residuals will be too small)" << endreq;
      m_UpdatorWarning = true;
      m_isUnbiased = 0;
    } else {
      m_isUnbiased = 0;
    }// end if no measured track parameter
  }
  return unbiasedTrkParameters;
}
**/
