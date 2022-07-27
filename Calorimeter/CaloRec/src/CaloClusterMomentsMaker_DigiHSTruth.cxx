/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/
//-----------------------------------------------------------------------
// File and Version Information:
//
// Description: see CaloClusterMomentsMaker.h
// 
// Environment:
//      Software developed for the ATLAS Detector at CERN LHC
//
// Author List:
//      Sven Menke
//
//-----------------------------------------------------------------------

#include "CaloClusterMomentsMaker_DigiHSTruth.h"
#include "CaloEvent/CaloClusterContainer.h"
#include "CaloEvent/CaloCluster.h"
#include "CaloGeoHelpers/proxim.h"
#include "CaloEvent/CaloPrefetch.h"
#include "CaloGeoHelpers/CaloPhiRange.h"
#include "CaloIdentifier/CaloCell_ID.h"
#include "AthAllocators/ArenaPoolSTLAllocator.h"

#include "GeoPrimitives/GeoPrimitives.h"
#include "GeoPrimitives/GeoPrimitivesHelpers.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CxxUtils/prefetch.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <sstream>


using CLHEP::deg;
using CLHEP::cm;

namespace {


  //FIXME, somehow make sure these names are in sync with the xAOD variable names
struct MomentName
{
  const char* name;
  xAOD::CaloCluster::MomentType mom;
};


// Must be sorted by name.
const MomentName moment_names[] = {
  { "ENERGY_DigiHSTruth",            xAOD::CaloCluster::ENERGY_DigiHSTruth },
  { "ETA_DigiHSTruth",               xAOD::CaloCluster::ETA_DigiHSTruth },
  { "PHI_DigiHSTruth",               xAOD::CaloCluster::PHI_DigiHSTruth },
  { "AVG_LAR_Q_DigiHSTruth",         xAOD::CaloCluster::AVG_LAR_Q_DigiHSTruth },
  { "AVG_TILE_Q_DigiHSTruth",        xAOD::CaloCluster::AVG_TILE_Q_DigiHSTruth },
  { "BADLARQ_FRAC_DigiHSTruth",      xAOD::CaloCluster::BADLARQ_FRAC_DigiHSTruth },
  { "BAD_CELLS_CORR_E_DigiHSTruth",  xAOD::CaloCluster::BAD_CELLS_CORR_E_DigiHSTruth },
  { "CELL_SIGNIFICANCE_DigiHSTruth", xAOD::CaloCluster::CELL_SIGNIFICANCE_DigiHSTruth },
  { "CELL_SIG_SAMPLING_DigiHSTruth", xAOD::CaloCluster::CELL_SIG_SAMPLING_DigiHSTruth },
  { "CENTER_LAMBDA_DigiHSTruth",     xAOD::CaloCluster::CENTER_LAMBDA_DigiHSTruth },
  { "CENTER_MAG_DigiHSTruth",        xAOD::CaloCluster::CENTER_MAG_DigiHSTruth },
  { "CENTER_X_DigiHSTruth",          xAOD::CaloCluster::CENTER_X_DigiHSTruth },
  { "CENTER_Y_DigiHSTruth",          xAOD::CaloCluster::CENTER_Y_DigiHSTruth },
  { "CENTER_Z_DigiHSTruth",          xAOD::CaloCluster::CENTER_Z_DigiHSTruth },
  { "DELTA_ALPHA_DigiHSTruth",       xAOD::CaloCluster::DELTA_ALPHA_DigiHSTruth },
  { "DELTA_PHI_DigiHSTruth",         xAOD::CaloCluster::DELTA_PHI_DigiHSTruth },
  { "DELTA_THETA_DigiHSTruth",       xAOD::CaloCluster::DELTA_THETA_DigiHSTruth },
  { "ENG_BAD_CELLS_DigiHSTruth",     xAOD::CaloCluster::ENG_BAD_CELLS_DigiHSTruth },
  { "ENG_BAD_HV_CELLS_DigiHSTruth",  xAOD::CaloCluster::ENG_BAD_HV_CELLS_DigiHSTruth },
  { "ENG_FRAC_CORE_DigiHSTruth",     xAOD::CaloCluster::ENG_FRAC_CORE_DigiHSTruth },
  { "ENG_FRAC_EM_DigiHSTruth",       xAOD::CaloCluster::ENG_FRAC_EM_DigiHSTruth },
  { "ENG_FRAC_MAX_DigiHSTruth",      xAOD::CaloCluster::ENG_FRAC_MAX_DigiHSTruth },
  { "ENG_POS_DigiHSTruth",           xAOD::CaloCluster::ENG_POS_DigiHSTruth },
  { "FIRST_ENG_DENS_DigiHSTruth",    xAOD::CaloCluster::FIRST_ENG_DENS_DigiHSTruth },
  { "FIRST_ETA_DigiHSTruth",         xAOD::CaloCluster::FIRST_ETA_DigiHSTruth },
  { "FIRST_PHI_DigiHSTruth",         xAOD::CaloCluster::FIRST_PHI_DigiHSTruth },
  { "ISOLATION_DigiHSTruth",         xAOD::CaloCluster::ISOLATION_DigiHSTruth },
  { "LATERAL_DigiHSTruth",           xAOD::CaloCluster::LATERAL_DigiHSTruth },
  { "LONGITUDINAL_DigiHSTruth",      xAOD::CaloCluster::LONGITUDINAL_DigiHSTruth },
  { "N_BAD_CELLS_DigiHSTruth",       xAOD::CaloCluster::N_BAD_CELLS_DigiHSTruth },
  { "N_BAD_HV_CELLS_DigiHSTruth",    xAOD::CaloCluster::N_BAD_HV_CELLS_DigiHSTruth },
  { "N_BAD_CELLS_CORR_DigiHSTruth",  xAOD::CaloCluster::N_BAD_CELLS_CORR_DigiHSTruth },
  { "SECOND_ENG_DENS_DigiHSTruth",   xAOD::CaloCluster::SECOND_ENG_DENS_DigiHSTruth },
  { "SECOND_LAMBDA_DigiHSTruth",     xAOD::CaloCluster::SECOND_LAMBDA_DigiHSTruth },
  { "SECOND_R_DigiHSTruth",          xAOD::CaloCluster::SECOND_R_DigiHSTruth },
  { "SIGNIFICANCE_DigiHSTruth",      xAOD::CaloCluster::SIGNIFICANCE_DigiHSTruth },
};

const MomentName* const moment_names_end =
  moment_names + sizeof(moment_names)/sizeof(moment_names[0]);

#if 0
bool operator< (const std::string& v, const MomentName& m)
{
  return strcmp (v.c_str(), m.name) < 0;
}
#endif
bool operator< (const MomentName& m, const std::string& v)
{
  return strcmp (m.name, v.c_str()) < 0;
}


}


//###############################################################################

CaloClusterMomentsMaker_DigiHSTruth::CaloClusterMomentsMaker_DigiHSTruth(const std::string& type, 
						 const std::string& name,
						 const IInterface* parent)
  :AthAlgTool(type, name, parent),
    m_calo_id(nullptr),
    m_maxAxisAngle(20*deg),
    m_minRLateral(4*cm),
    m_minLLongitudinal(10*cm),
    m_minBadLArQuality(4000),
    m_calculateSignificance(false),
    m_calculateIsolation(false),
    m_calculateLArHVFraction(false),
    m_twoGaussianNoise(false),
    m_caloDepthTool("CaloDepthTool",this),
    m_larHVFraction("LArHVFraction",this),
    m_absOpt(false) 
{
  // Name(s) of Moments to calculate
  declareProperty("MomentsNames",m_momentsNames);

  // Name(s) of Moments which can be stored on the AOD - all others go to ESD
  // m_momentsNamesAOD.push_back(std::string("FIRST_PHI"));
  // m_momentsNamesAOD.push_back(std::string("FIRST_ETA"));
  // m_momentsNamesAOD.push_back(std::string("SECOND_R"));
  // m_momentsNamesAOD.push_back(std::string("SECOND_LAMBDA"));
  // m_momentsNamesAOD.push_back(std::string("CENTER_LAMBDA"));
  // m_momentsNamesAOD.push_back(std::string("FIRST_ENG_DENS"));
  // m_momentsNamesAOD.push_back(std::string("ENG_BAD_CELLS"));
  // m_momentsNamesAOD.push_back(std::string("N_BAD_CELLS"));

  //declareProperty("AODMomentsNames",m_momentsNamesAOD);
  // maximum allowed angle between shower axis and the vector pointing
  // to the shower center from the IP in degrees. This property is need
  // to protect against cases where all significant cells are in one sampling
  // and the shower axis can not be defined from them
  declareProperty("MaxAxisAngle",m_maxAxisAngle);
  declareProperty("MinRLateral",m_minRLateral);
  declareProperty("MinLLongitudinal",m_minLLongitudinal);
  declareProperty("MinBadLArQuality",m_minBadLArQuality);
  // use 2-gaussian noise for Tile
  declareProperty("TwoGaussianNoise",m_twoGaussianNoise);
  declareProperty("LArHVFraction",m_larHVFraction,"Tool Handle for LArHVFraction");

  /// Not used anymore (with xAOD), but required when configured from COOL.
  declareProperty("AODMomentsNames",m_momentsNamesAOD);
  // Use weighting of neg. clusters option?
  declareProperty("WeightingOfNegClusters", m_absOpt);
}

//###############################################################################

StatusCode CaloClusterMomentsMaker_DigiHSTruth::initialize()
{

  //FIXME: All that could be done at initialize!
  m_calculateSignificance = false;
  m_calculateIsolation = false;

  // translate all moment names specified in MomentsNames property to moment enums,
  // check that they are all valid and there are no repeating names
  for(const auto& name: m_momentsNames) {
    const MomentName* it =
      std::lower_bound (moment_names, moment_names_end, name);
      //std::find(moment_names, moment_names_end, name);
      //moment_names.find(name);

		for(const auto& testName: moment_names){
			if(name == testName.name) it = &testName;
		}
    if (it != moment_names_end) {
      m_validMoments.push_back (it->mom);
      switch (it->mom) {
      case xAOD::CaloCluster::SIGNIFICANCE_DigiHSTruth:
      case xAOD::CaloCluster::CELL_SIGNIFICANCE_DigiHSTruth:
        m_calculateSignificance = true;
        break;
      case xAOD::CaloCluster::ISOLATION_DigiHSTruth:
        m_calculateIsolation = true;
        break;
      case xAOD::CaloCluster::ENG_BAD_HV_CELLS_DigiHSTruth:
        m_calculateLArHVFraction = true;
      default:
        break;
      }
    }
    else {
      msg(MSG::ERROR) << "Moment " << name
		      << " is not a valid Moment name and will be ignored! "
		      << "Valid names are:";
      int count = 0;
      for (const MomentName& m : moment_names)
	msg() << ((count++)==0?" ":", ") << m.name;
      msg() << endmsg;
    }
  }

  // sort and remove duplicates, order is not required for any of the code below
  // but still may be useful property
  std::stable_sort(m_validMoments.begin(), m_validMoments.end());
/*
  m_validMoments.erase(std::unique(m_validMoments.begin(),
                                   m_validMoments.end()),
                       m_validMoments.end());
*/

  
  /*
  // translate moment names in AODMomentsNames property into set of enums,
  // only take valid names which are also in MomentsNames property
  m_momentsAOD.reserve(m_momentsNamesAOD.size());
  for(const auto& name: m_momentsNamesAOD) {
    const MomentName* it =
      std::lower_bound (moment_names, moment_names_end, name);
    if (it != moment_names_end) {
      if (std::find(m_validMoments.begin(), m_validMoments.end(), it->mom)
          != m_validMoments.end())
      {
        m_momentsAOD.push_back(it->mom);
      }
    }
  }
  */
  
  ATH_CHECK(detStore()->retrieve(m_calo_id,"CaloCell_ID"));
  
  ATH_CHECK(m_caloDepthTool.retrieve());
  ATH_CHECK(m_caloMgrKey.initialize());

  if (m_calculateSignificance) {
    ATH_CHECK(m_noiseCDOKey.initialize());
  }

  if (m_calculateLArHVFraction) {
    ATH_CHECK(m_larHVFraction.retrieve());
  }
  else {
    m_larHVFraction.disable();
  }

  return StatusCode::SUCCESS;
  
}

//#############################################################################

namespace CaloClusterMomentsMaker_detail {

struct cellinfo {
  double x;
  double y;
  double z;
  double energy;
  double eta;
  double phi;
  double r;
  double lambda;
  double volume;
  CaloCell_ID::CaloSample sample;
};

} // namespace CaloClusterMomentsMaker_detail

StatusCode
CaloClusterMomentsMaker_DigiHSTruth::execute(const EventContext& ctx,
                                 xAOD::CaloClusterContainer *theClusColl)
  const
{

  ATH_MSG_DEBUG("Executing " << name());

  SG::ReadHandle<CaloCellContainer> signalCells(m_signalCellKey,ctx);

  // Maps cell IdentifierHash to cluster index in cluster collection.
  // Only used when cluster isolation moment is calculated.
  using clusterIdx_t = std::uint16_t;
  typedef std::pair<clusterIdx_t, clusterIdx_t> clusterPair_t;
  std::vector<clusterPair_t> clusterIdx;
  const clusterIdx_t noCluster = std::numeric_limits<clusterIdx_t>::max();

  const CaloNoise* noise=nullptr;
  if (m_calculateSignificance) {
    SG::ReadCondHandle<CaloNoise> noiseHdl{m_noiseCDOKey,ctx};
    noise=*noiseHdl;
  }

  SG::ReadCondHandle<CaloDetDescrManager> caloMgrHandle{ m_caloMgrKey, ctx };
  const CaloDetDescrManager* caloDDMgr = *caloMgrHandle;
 
  // Counters for number of empty and non-empty neighbor cells per sampling layer
  // Only used when cluster isolation moment is calculated.
  int nbEmpty[CaloCell_ID::Unknown];
  int nbNonEmpty[CaloCell_ID::Unknown];

  // prepare stuff from entire collection in case isolation moment
  // should be calculated

  if ( m_calculateIsolation ) {

    if (theClusColl->size() >= noCluster) {
      msg(MSG::ERROR) << "Too many clusters" << endmsg;
      return StatusCode::FAILURE;
    }

    // initialize with "empty" values
    clusterIdx.resize(m_calo_id->calo_cell_hash_max(),
                      clusterPair_t(noCluster, noCluster));

    int iCluster = 0;
    for (xAOD::CaloCluster* theCluster : *theClusColl) {
      // loop over all cell members and fill cell vector for used cells
      xAOD::CaloCluster::cell_iterator cellIter    = theCluster->cell_begin();
      xAOD::CaloCluster::cell_iterator cellIterEnd = theCluster->cell_end();
      for(; cellIter != cellIterEnd; cellIter++ ){
        CxxUtils::prefetchNext(cellIter, cellIterEnd);
        const CaloCell* myCell = *cellIter;

        const CaloDetDescrElement * caloDDE = myCell->caloDDE();
        IdentifierHash hashid=caloDDE->calo_hash() ;
        if(! hashid.is_valid() ) continue;
        if(hashid >= (signalCells)->size()) continue;
        myCell = (*signalCells).findCell(hashid);
        if(!myCell) continue;


	Identifier myId = myCell->ID();
	IdentifierHash myHashId = m_calo_id->calo_cell_hash(myId); 
	if ( clusterIdx[(unsigned int)myHashId].first != noCluster) {
	  // check weight and assign to current cluster if weight is > 0.5
	  double weight = cellIter.weight();
	  if ( weight > 0.5 )
	    clusterIdx[(unsigned int)myHashId].first = iCluster;
	}
	else {
	  clusterIdx[(unsigned int)myHashId].first = iCluster;
	}
      }
      ++iCluster;
    }
  }

  // Move allocation of temporary arrays outside the cluster loop.
  // That way, we don't need to delete and reallocate them
  // each time through the loop.

  std::vector<CaloClusterMomentsMaker_detail::cellinfo> cellinfo;
  std::vector<double> maxSampE (CaloCell_ID::Unknown);
  std::vector<double> myMoments(m_validMoments.size(),0);
  std::vector<double> myNorms(m_validMoments.size(),0);
  std::vector<IdentifierHash> theNeighbors;
  // loop over individual clusters
  xAOD::CaloClusterContainer::iterator clusIter = theClusColl->begin();
  xAOD::CaloClusterContainer::iterator clusIterEnd = theClusColl->end();
  int iClus = 0;
  for( ;clusIter!=clusIterEnd;++clusIter,++iClus) {
    xAOD::CaloCluster * theCluster = *clusIter;

    double w(0),xc(0),yc(0),zc(0);
    double eBad(0),ebad_dac(0),ePos(0),eBadLArQ(0),sumSig2(0),maxAbsSig(0);
    double eLAr2(0),eLAr2Q(0);
    double eTile2(0),eTile2Q(0);
    double eBadLArHV(0);
    int nbad(0),nbad_dac(0),nBadLArHV(0);
    unsigned int ncell(0),i,nSigSampl(0);
    unsigned int theNumOfCells = theCluster->size();
		double theClusterEnergy = 0;
		double theClusterAbsEnergy = 0;
		double theClusterEta = 0;
		double theClusterPhi = 0;
    
    // these two are needed for the LATERAL moment
    int iCellMax(-1);
    int iCellScndMax(-1);

    if (cellinfo.capacity() == 0)
      cellinfo.reserve (theNumOfCells*2);
    cellinfo.resize (theNumOfCells);
    
    double phi0 = theCluster->phi();;
    for(i=0;i<(unsigned int)CaloCell_ID::Unknown;i++) 
      maxSampE[i] = 0;
    
    if ( !m_momentsNames.empty() ) {
      std::fill (myMoments.begin(), myMoments.end(), 0);
      std::fill (myNorms.begin(),   myNorms.end(),   0);
      if ( m_calculateIsolation ) {
        std::fill_n(nbNonEmpty, CaloCell_ID::Unknown, 0);
        std::fill_n(nbEmpty, CaloCell_ID::Unknown, 0);
      }
      
      // loop over all cell members and calculate the center of mass
      xAOD::CaloCluster::cell_iterator cellIter    = theCluster->cell_begin();
      xAOD::CaloCluster::cell_iterator cellIterEnd = theCluster->cell_end();
      for(; cellIter != cellIterEnd; cellIter++ ){
        CaloPrefetch::nextDDE(cellIter, cellIterEnd);
        const CaloCell* pCell = (*cellIter);
        const CaloDetDescrElement * caloDDE = pCell->caloDDE();
        IdentifierHash hashid=caloDDE->calo_hash() ;
        if(! hashid.is_valid() ) continue;

        if(hashid >= (signalCells)->size()) continue;
        pCell = (*signalCells).findCell(hashid);

      	Identifier myId = pCell->ID();

      	const CaloDetDescrElement* myCDDE = pCell->caloDDE();
      	double ene = pCell->e();
        if(m_absOpt) ene = std::abs(ene);
      	double weight = cellIter.weight();//theCluster->getCellWeight(cellIter);

        double thePhi;
        double cellPhi = myCDDE->phi();
        thePhi = proxim (cellPhi, phi0);

        theClusterEnergy += weight * ene;
        theClusterAbsEnergy += weight*std::abs(ene);
        theClusterEta += weight*std::abs(ene)*pCell->eta(); 
        theClusterPhi += weight*std::abs(ene)* thePhi; 
      	if ( pCell->badcell() )       {
      	  eBad += ene*weight;
      	  nbad++;
      	  if(ene!=0){
      	    ebad_dac+=ene*weight;
      	    nbad_dac++;
	        }
      	}
  	else {
	  if ( ! (myCDDE->is_tile()) 
	       && ((pCell->provenance() & 0x2000) == 0x2000) 
	       && !((pCell->provenance() & 0x0800) == 0x0800)) {
	    if ( pCell->quality() > m_minBadLArQuality ) {
	      eBadLArQ += ene*weight;
	    }
	    eLAr2  += ene*weight*ene*weight;
	    eLAr2Q += ene*weight*ene*weight*pCell->quality();
	  }
	  if ( myCDDE->is_tile() ) {
	    uint16_t tq = pCell->quality();
	    uint8_t tq1 = (0xFF00&tq)>>8; // quality in channel 1
	    uint8_t tq2 = (0xFF&tq); // quality in channel 2
	    // reject cells with either 0xFF00 or 0xFF
	    if ( ((tq1&0xFF) != 0xFF) && ((tq2&0xFF) != 0xFF) ) {    
	      eTile2  += ene*weight*ene*weight;
	      // take the worse of both qualities (one might be 0 in
	      // 1-channel cases)
	      eTile2Q += ene*weight*ene*weight*(tq1>tq2?tq1:tq2);
	    }
	  } 
	}
	if ( ene > 0 ) {
	  ePos += ene*weight;
	}
	if ( m_calculateSignificance ) {
	  const float sigma = m_twoGaussianNoise ?			\
	    noise->getEffectiveSigma(pCell->ID(),pCell->gain(),pCell->energy()) : \
	    noise->getNoise(pCell->ID(),pCell->gain());
	  sumSig2 += sigma*sigma;
	  // use geomtery weighted energy of cell for leading cell significance
	  double Sig = (sigma>0?ene*weight/sigma:0);
	  if ( std::abs(Sig) > std::abs(maxAbsSig) ) {
	    maxAbsSig = Sig;
            nSigSampl = myCDDE->getSampling();
	  }
	}
	if ( m_calculateIsolation ) {
	  // get all 2D Neighbours if the cell is not inside another cluster with
	  // larger weight

	  IdentifierHash myHashId = m_calo_id->calo_cell_hash(myId);
	  if ( clusterIdx[myHashId].first == iClus ) {
            theNeighbors.clear();
	    m_calo_id->get_neighbours(myHashId, LArNeighbours::all2D, theNeighbors);
	    for (const auto& nhash: theNeighbors) {
              clusterPair_t& idx = clusterIdx[nhash];

	      // only need to look at each cell once per cluster
	      if ( idx.second == iClus ) continue;
	      idx.second = iClus;

	      if ( idx.first == noCluster ) {
	        ++ nbEmpty[m_calo_id->calo_sample(nhash)];
	      } else if ( idx.first != iClus ) {
                ++ nbNonEmpty[m_calo_id->calo_sample(nhash)];
	      }

	    }
	  }
	}

	if ( myCDDE && ene > 0. && weight > 0) {
	  // get all geometric information needed ...
          CaloClusterMomentsMaker_detail::cellinfo& ci = cellinfo[ncell];
	  ci.x      = myCDDE->x();
	  ci.y      = myCDDE->y();
	  ci.z      = myCDDE->z();
	  ci.eta    = myCDDE->eta();
	  ci.phi    = myCDDE->phi();
	  ci.energy = ene*weight;
	  ci.volume = myCDDE->volume();
	  ci.sample = myCDDE->getSampling();
	  if ( ci.energy > maxSampE[(unsigned int)ci.sample] )
	    maxSampE[(unsigned int)ci.sample] = ci.energy;
	  
	  if (iCellMax < 0 || ci.energy > cellinfo[iCellMax].energy ) {
	    iCellScndMax = iCellMax;
	    iCellMax = ncell;
	  }
	  else if (iCellScndMax < 0 ||
                   ci.energy > cellinfo[iCellScndMax].energy )
          {
	    iCellScndMax = ncell;
	  }
	  
	  xc += ci.energy*ci.x;
	  yc += ci.energy*ci.y;
	  zc += ci.energy*ci.z;
	  w  += ci.energy;
	 
	  ncell++;
	}
      } //end of loop over all cells

      if ( w > 0 ) {
	xc/=w;
	yc/=w;
	zc/=w;
  Amg::Vector3D showerCenter(xc,yc,zc);
	w=0;
	

	//log << MSG::WARNING << "Found bad cells " <<  xbad_dac << " " << ybad_dac << " " << zbad_dac << " " << ebad_dac <<  endmsg;
	//log << MSG::WARNING << "Found Cluster   " <<  xbad_dac << " " << ybad_dac << " " << zbad_dac << " " <<  endmsg;
	// shower axis is just the vector pointing from the IP to the shower center
	// in case there are less than 3 cells in the cluster
	
  Amg::Vector3D showerAxis(xc,yc,zc);
  Amg::setMag(showerAxis,1.0);

	// otherwise the principal direction with the largest absolute 
	// eigenvalue will be used unless it's angle w.r.t. the vector pointing
	// from the IP to the shower center is larger than allowed by the 
	// property m_maxAxisAngle

	double angle(0),deltaPhi(0),deltaTheta(0);
	if ( ncell > 2 ) {
	  Eigen::Matrix3d C=Eigen::Matrix3d::Zero();
	  for(i=0;i<ncell;i++) {
            const CaloClusterMomentsMaker_detail::cellinfo& ci = cellinfo[i];
            const double e2 = ci.energy * ci.energy;
            
	    C(0,0) += e2*(ci.x-xc)*(ci.x-xc);
	    C(1,0) += e2*(ci.x-xc)*(ci.y-yc);
	    C(2,0) += e2*(ci.x-xc)*(ci.z-zc);
	    
	    C(1,1) += e2*(ci.y-yc)*(ci.y-yc);
	    C(2,1) += e2*(ci.y-yc)*(ci.z-zc);
	    
	    C(2,2) += e2*(ci.z-zc)*(ci.z-zc);
	    w += e2;
	  } 

	  C/=w;
	  

	  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(C);
	  if (eigensolver.info() != Eigen::Success) {
	    msg(MSG::WARNING) << "Failed to compute Eigenvalues -> Can't determine shower axis" << endmsg;
	  }
	  else {
	    // don't use the principal axes if at least one of the 3 
	    // diagonal elements is 0


      const Eigen::Vector3d& S=eigensolver.eigenvalues();
      const Eigen::Matrix3d& U=eigensolver.eigenvectors();

            const double epsilon = 1.E-6;

      if ( !( std::abs(S[0]) < epsilon || std::abs(S[1]) < epsilon || std::abs(S[2]) < epsilon ) ) {

        Amg::Vector3D prAxis(showerAxis);

	      int iEigen = -1;
	    
	      for (i=0;i<3;i++) {
      Amg::Vector3D tmpAxis=U.col(i);

      // calculate the angle    
      double tmpAngle=Amg::angle(tmpAxis,showerAxis);

      if ( tmpAngle > 90*deg ) {
        tmpAngle = 180*deg - tmpAngle;
        tmpAxis = -tmpAxis;
      }

      if ( iEigen == -1 || tmpAngle < angle ) {
        iEigen = i;
        angle = tmpAngle;
        prAxis = tmpAxis;
      }

	      }//end for loop 	  
	    
	      // calculate theta and phi angle differences
	      
	      deltaPhi = CaloPhiRange::diff(showerAxis.phi(),prAxis.phi());
	    
	      deltaTheta = showerAxis.theta() - prAxis.theta();
	    
	      // check the angle
	    
	      if ( angle < m_maxAxisAngle ) {
		showerAxis = prAxis;
        }
	      else 
    ATH_MSG_DEBUG("principal Direction (" << prAxis[Amg::x] << ", "
            << prAxis[Amg::y] << ", " << prAxis[Amg::z] << ") deviates more than "
            << m_maxAxisAngle*(1./deg)
            << " deg from IP-to-ClusterCenter-axis (" << showerAxis[Amg::x] << ", "
            << showerAxis[Amg::y] << ", " << showerAxis[Amg::z] << ")");
	    }//end if !S[i]==0
	  }// end else got Eigenvalues
	} //end if ncell>2
      
	ATH_MSG_DEBUG("Shower Axis = (" << showerAxis.x() << ", "
		      << showerAxis.y() << ", " << showerAxis.z() << ")");
	
	// calculate radial distance from and the longitudinal distance
	// along the shower axis for each cell. The cluster center is 
	// at r=0 and lambda=0
	
  for (auto& ci : cellinfo) {
    const Amg::Vector3D currentCell(ci.x,ci.y,ci.z);
	  // calculate distance from shower axis r
	  ci.r = ((currentCell-showerCenter).cross(showerAxis)).mag();
	  // calculate distance from shower center along shower axis
	  ci.lambda = (currentCell-showerCenter).dot(showerAxis);
	}  
	
	// loop over all positive energy cells and calculate all desired moments
	
	// define common norm for all simple moments
	double commonNorm = 0;
        double phi0 = ncell > 0 ? cellinfo[0].phi : 0;
	for(i=0;i<ncell;i++) {
          const CaloClusterMomentsMaker_detail::cellinfo& ci = cellinfo[i];
	  // loop over all valid moments
	  commonNorm += ci.energy;
	  for(size_t iMoment = 0, size = m_validMoments.size();
              iMoment != size;
              ++ iMoment)
          {
	    // now calculate the actual moments
	    switch (m_validMoments[iMoment]) {
       case xAOD::CaloCluster::FIRST_ETA_DigiHSTruth:
        myMoments[iMoment] += ci.energy*ci.eta;
        break;
	    case xAOD::CaloCluster::FIRST_PHI_DigiHSTruth:
	      // first cell decides the sign in order to avoid
	      // overlap problem at phi = -pi == +pi
	      // need to be normalized to the range [-pi,+pi] in the end
        myMoments[iMoment] += ci.energy * proxim (ci.phi, phi0);
	      break;
	    case xAOD::CaloCluster::SECOND_R_DigiHSTruth:
	      myMoments[iMoment] += ci.energy*ci.r*ci.r;
	      break;
	    case xAOD::CaloCluster::SECOND_LAMBDA_DigiHSTruth:
	      myMoments[iMoment] += ci.energy*ci.lambda*ci.lambda;
	      break;
	    case xAOD::CaloCluster::LATERAL_DigiHSTruth:
	      if ( (int)i != iCellMax && (int)i != iCellScndMax ) {
		myMoments[iMoment] += ci.energy*ci.r*ci.r;
		myNorms[iMoment] += ci.energy*ci.r*ci.r;
	      }
	      else {
		double rm = ci.r;
		if ( rm < m_minRLateral ) 
		  rm = m_minRLateral;
		myNorms[iMoment] += rm*rm*ci.energy;
	      }
	      break;
	    case xAOD::CaloCluster::LONGITUDINAL_DigiHSTruth:
	      if ( (int)i != iCellMax && (int)i != iCellScndMax ) {
		myMoments[iMoment] += ci.energy*ci.lambda*ci.lambda;
		myNorms[iMoment] += ci.energy*ci.lambda*ci.lambda;
	      }
	      else {
		double lm = ci.lambda;
		if ( lm < m_minLLongitudinal ) 
		  lm = m_minLLongitudinal;
		myNorms[iMoment] += lm*lm*ci.energy;
	      }
	      break;
	    case xAOD::CaloCluster::FIRST_ENG_DENS_DigiHSTruth:
	      if ( ci.volume > 0 ) {
		myMoments[iMoment] += ci.energy*ci.energy/ci.volume;
		myNorms[iMoment] += ci.energy;
	      }
	      break;
	    case xAOD::CaloCluster::SECOND_ENG_DENS_DigiHSTruth:
	      if ( ci.volume > 0 ) {
		myMoments[iMoment] += ci.energy*std::pow(ci.energy/ci.volume,2);
		myNorms[iMoment] += ci.energy;
	      }
	      break;
	    case xAOD::CaloCluster::ENG_FRAC_EM_DigiHSTruth:
	      if ( ci.sample == CaloCell_ID::EMB1 
		   || ci.sample == CaloCell_ID::EMB2 
		   || ci.sample == CaloCell_ID::EMB3 
		   || ci.sample == CaloCell_ID::EME1
		   || ci.sample == CaloCell_ID::EME2
		   || ci.sample == CaloCell_ID::EME3
		   || ci.sample == CaloCell_ID::FCAL0 )
		myMoments[iMoment] += ci.energy;
	      break;
	    case xAOD::CaloCluster::ENG_FRAC_MAX_DigiHSTruth:
	      if ( (int)i == iCellMax ) 
		myMoments[iMoment] = ci.energy;
	      break;
	    default:
	      // nothing to be done for other moments
	      break;
	    }
	  }
	} //end of loop over cell

	const auto hvFrac=m_larHVFraction->getLArHVFrac(theCluster->getCellLinks(),ctx);
	eBadLArHV= hvFrac.first;
	nBadLArHV=hvFrac.second;

	
	// assign moments which don't need the loop over the cells
        for (size_t iMoment = 0, size = m_validMoments.size();
             iMoment != size;
             ++ iMoment)
        {
	  // now calculate the actual moments
          switch (m_validMoments[iMoment]) {
	  case xAOD::CaloCluster::FIRST_ETA_DigiHSTruth:
	  case xAOD::CaloCluster::FIRST_PHI_DigiHSTruth:
	  case xAOD::CaloCluster::SECOND_R_DigiHSTruth:
	  case xAOD::CaloCluster::SECOND_LAMBDA_DigiHSTruth:
	  case xAOD::CaloCluster::ENG_FRAC_EM_DigiHSTruth:
	  case xAOD::CaloCluster::ENG_FRAC_MAX_DigiHSTruth:
	    myNorms[iMoment] = commonNorm;
	    break;
	  case xAOD::CaloCluster::DELTA_PHI_DigiHSTruth:
	    myMoments[iMoment] = deltaPhi;
	    break;
	  case xAOD::CaloCluster::DELTA_THETA_DigiHSTruth:
	    myMoments[iMoment] = deltaTheta;
	    break;
	  case xAOD::CaloCluster::DELTA_ALPHA_DigiHSTruth:
	    myMoments[iMoment] = angle;
	    break;
	  case xAOD::CaloCluster::CENTER_X_DigiHSTruth:
	    myMoments[iMoment] = showerCenter.x();
	    break;
	  case xAOD::CaloCluster::CENTER_Y_DigiHSTruth:
	    myMoments[iMoment] = showerCenter.y();
	    break;
	  case xAOD::CaloCluster::CENTER_Z_DigiHSTruth:
	    myMoments[iMoment] = showerCenter.z();
	    break;
	  case xAOD::CaloCluster::CENTER_MAG_DigiHSTruth:
	    myMoments[iMoment] = showerCenter.mag();
	    break;
	  case xAOD::CaloCluster::CENTER_LAMBDA_DigiHSTruth:
	    // calculate the longitudinal distance along the shower axis
	    // of the shower center from the calorimeter start
	    
	    // first need calo boundary at given eta phi try LAREM barrel
	    // first, then LAREM endcap OW, then LAREM endcap IW, then
	    // FCal
	    {
	      double r_calo(0),z_calo(0),lambda_c(0); 
	      r_calo = m_caloDepthTool->get_entrance_radius(CaloCell_ID::EMB1,
							    showerCenter.eta(),
							    showerCenter.phi(),
                  caloDDMgr);
	      if ( r_calo == 0 ) {
		z_calo = m_caloDepthTool->get_entrance_z(CaloCell_ID::EME1,
							 showerCenter.eta(),
							 showerCenter.phi(),
               caloDDMgr);
		if ( z_calo == 0 ) 
		  z_calo = m_caloDepthTool->get_entrance_z(CaloCell_ID::EME2,
							   showerCenter.eta(),
							   showerCenter.phi(),
                 caloDDMgr);
		if ( z_calo == 0 ) 
		  z_calo = m_caloDepthTool->get_entrance_z(CaloCell_ID::FCAL0,
							   showerCenter.eta(),
							   showerCenter.phi(),
                 caloDDMgr);
		if ( z_calo == 0 ) // for H6 TB without EMEC outer wheel 
		  z_calo = m_caloDepthTool->get_entrance_z(CaloCell_ID::HEC0,
							   showerCenter.eta(),
							   showerCenter.phi(),
                 caloDDMgr);
		if ( z_calo != 0 && showerAxis.z() != 0 ) {
		  lambda_c = std::abs((z_calo-showerCenter.z())/showerAxis.z());
		}
	      }
	      else {
		double r_s2 = showerAxis.x()*showerAxis.x()
		  +showerAxis.y()*showerAxis.y();
		double r_cs = showerAxis.x()*showerCenter.x()
		  +showerAxis.y()*showerCenter.y();
		double r_cr = showerCenter.x()*showerCenter.x()
		  +showerCenter.y()*showerCenter.y()-r_calo*r_calo;
		if ( r_s2 > 0 ) {
		  double det = r_cs*r_cs/(r_s2*r_s2) - r_cr/r_s2;
		  if ( det > 0 ) {
		    det = sqrt(det);
		    double l1(-r_cs/r_s2);
		    double l2(l1);
		    l1 += det;
		    l2 -= det;
		    if ( std::abs(l1) < std::abs(l2) ) 
		      lambda_c = std::abs(l1);
		    else
		      lambda_c = std::abs(l2);
		  }
		}
	      }
	      myMoments[iMoment] = lambda_c;
	    }
	    break;
	  case xAOD::CaloCluster::ENG_FRAC_CORE_DigiHSTruth:
	    for(i=0;i<(int)CaloCell_ID::Unknown;i++) 
	      myMoments[iMoment] += maxSampE[i];
	    myNorms[iMoment] = commonNorm;
	    break;
	  case xAOD::CaloCluster::ISOLATION_DigiHSTruth:
	    {
	      // loop over empty and filled perimeter cells and
	      // get a weighted ratio by means of energy fraction per layer
	      for(unsigned int i=0; i != CaloSampling::Unknown; ++ i) {
		const xAOD::CaloCluster::CaloSample s=( xAOD::CaloCluster::CaloSample)(i);
		if (theCluster->hasSampling(s)) {
		  const double eSample = theCluster->eSample(s);
		  if (eSample > 0) {
		    int nAll = nbEmpty[i]+nbNonEmpty[i];
		    if (nAll > 0) {
		      myMoments[iMoment] += (eSample*nbEmpty[i])/nAll;
		      myNorms[iMoment] += eSample;
		    }
		  }//end of eSample>0
		}//end has sampling
	      }//end loop over samplings
	    }
	    break;
	  case xAOD::CaloCluster::ENG_BAD_CELLS_DigiHSTruth:
	    myMoments[iMoment] = eBad;
	    break;
	  case xAOD::CaloCluster::N_BAD_CELLS_DigiHSTruth:
	    myMoments[iMoment] = nbad;
            break;
          case xAOD::CaloCluster::N_BAD_CELLS_CORR_DigiHSTruth:
            myMoments[iMoment] = nbad_dac;
            break;
	  case xAOD::CaloCluster::BAD_CELLS_CORR_E_DigiHSTruth:
	    myMoments[iMoment] = ebad_dac;
            break;
	  case xAOD::CaloCluster::BADLARQ_FRAC_DigiHSTruth:
	    myMoments[iMoment] = eBadLArQ/(theCluster->e()!=0.?theCluster->e():1.);
            break;
	  case xAOD::CaloCluster::ENG_POS_DigiHSTruth:
	    myMoments[iMoment] = ePos;
            break;
	  case xAOD::CaloCluster::SIGNIFICANCE_DigiHSTruth:
	    myMoments[iMoment] = (sumSig2>0?theCluster->e()/sqrt(sumSig2):0.);
            break;
	  case xAOD::CaloCluster::CELL_SIGNIFICANCE_DigiHSTruth:
	    myMoments[iMoment] = maxAbsSig;
            break;
	  case xAOD::CaloCluster::CELL_SIG_SAMPLING_DigiHSTruth:
	    myMoments[iMoment] = nSigSampl;
            break;
	  case xAOD::CaloCluster::AVG_LAR_Q_DigiHSTruth:
	    myMoments[iMoment] = eLAr2Q/(eLAr2>0?eLAr2:1);
            break;
	  case xAOD::CaloCluster::AVG_TILE_Q_DigiHSTruth:
	    myMoments[iMoment] = eTile2Q/(eTile2>0?eTile2:1);
            break;
	  case xAOD::CaloCluster::ENG_BAD_HV_CELLS_DigiHSTruth:
	    myMoments[iMoment] = eBadLArHV;
	    break;
	  case xAOD::CaloCluster::N_BAD_HV_CELLS_DigiHSTruth:
	    myMoments[iMoment] = nBadLArHV;
            break;
	  case xAOD::CaloCluster::ENERGY_DigiHSTruth:
	    myMoments[iMoment] = theClusterEnergy;
            break;
	  case xAOD::CaloCluster::ETA_DigiHSTruth:
      if(theClusterAbsEnergy > 0)
	    myMoments[iMoment] = theClusterEta / theClusterAbsEnergy;
      else{
	    myMoments[iMoment] = 0;
			}
            break;
	  case xAOD::CaloCluster::PHI_DigiHSTruth:
      if(theClusterAbsEnergy > 0)
	    myMoments[iMoment] = CaloPhiRange::fix(theClusterPhi / theClusterAbsEnergy);
      else{
	    myMoments[iMoment] = 0;
			}
            break;
	  default:
	    // nothing to be done for other moments
	    break;
	  }
	}
      }
      
      // normalize moments and copy to Cluster Moment Store
      size_t size= m_validMoments.size();
      for (size_t iMoment = 0; iMoment != size; ++iMoment) {
        xAOD::CaloCluster::MomentType moment = m_validMoments[iMoment];
	if ( myNorms[iMoment] != 0 ) 
	  myMoments[iMoment] /= myNorms[iMoment];
	if ( moment == xAOD::CaloCluster::FIRST_PHI_DigiHSTruth ) 
	  myMoments[iMoment] = CaloPhiRange::fix(myMoments[iMoment]);
	
	theCluster->insertMoment(moment,myMoments[iMoment]);
      }
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloClusterMomentsMaker_DigiHSTruth::finalize()
{
  return StatusCode::SUCCESS;
}
