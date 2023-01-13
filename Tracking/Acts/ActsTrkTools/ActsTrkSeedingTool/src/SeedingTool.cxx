//Dear emacs, this is -*- c++ -*-
/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "src/SeedingTool.h"

// ACTS
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"

namespace ActsTrk {
   SeedingTool::SeedingTool(const std::string& type,
    const std::string& name,
    const IInterface* parent)
    : base_class(type, name, parent)
  {}
  
  StatusCode SeedingTool::initialize() {
    ATH_MSG_DEBUG("Initializing " << name() << "...");

    ATH_MSG_DEBUG("Properties Summary:");
    ATH_MSG_DEBUG("   " << m_zBinNeighborsTop);
    ATH_MSG_DEBUG("   " << m_zBinNeighborsBottom);
    ATH_MSG_DEBUG("   " << m_numPhiNeighbors);

    ATH_MSG_DEBUG(" *  Used by SpacePointGridConfig");
    ATH_MSG_DEBUG("   " << m_minPt);
    ATH_MSG_DEBUG("   " << m_cotThetaMax);
    ATH_MSG_DEBUG("   " << m_impactMax);
    ATH_MSG_DEBUG("   " << m_zMin);
    ATH_MSG_DEBUG("   " << m_zMax);
    ATH_MSG_DEBUG("   " << m_gridPhiMin);
    ATH_MSG_DEBUG("   " << m_gridPhiMax);
    ATH_MSG_DEBUG("   " << m_zBinEdges);
    ATH_MSG_DEBUG("   " << m_deltaRMax);
    ATH_MSG_DEBUG("   " << m_gridRMax);
    ATH_MSG_DEBUG("   " << m_bFieldInZ);
    ATH_MSG_DEBUG("   " << m_phiBinDeflectionCoverage);

    ATH_MSG_DEBUG(" * Used by SeedFinderConfig:");
    ATH_MSG_DEBUG("   " << m_minPt);
    ATH_MSG_DEBUG("   " << m_cotThetaMax);
    ATH_MSG_DEBUG("   " << m_impactMax);
    ATH_MSG_DEBUG("   " << m_zMin);
    ATH_MSG_DEBUG("   " << m_zMax);
    ATH_MSG_DEBUG("   " << m_zBinEdges);
    ATH_MSG_DEBUG("   " << m_rMax);
    ATH_MSG_DEBUG("   " << m_deltaRMin);
    ATH_MSG_DEBUG("   " << m_deltaRMax);
    ATH_MSG_DEBUG("   " << m_deltaRMinTopSP);
    ATH_MSG_DEBUG("   " << m_deltaRMaxTopSP);
    ATH_MSG_DEBUG("   " << m_deltaRMinBottomSP);
    ATH_MSG_DEBUG("   " << m_deltaRMaxBottomSP);
    ATH_MSG_DEBUG("   " << m_deltaZMax);
    ATH_MSG_DEBUG("   " << m_collisionRegionMin);
    ATH_MSG_DEBUG("   " << m_collisionRegionMax);
    ATH_MSG_DEBUG("   " << m_sigmaScattering);
    ATH_MSG_DEBUG("   " << m_maxPtScattering);
    ATH_MSG_DEBUG("   " << m_radLengthPerSeed);
    ATH_MSG_DEBUG("   " << m_maxSeedsPerSpM);
    ATH_MSG_DEBUG("   " << m_interactionPointCut);
    ATH_MSG_DEBUG("   " << m_arithmeticAverageCotTheta);
    ATH_MSG_DEBUG("   " << m_skipPreviousTopSP);
    ATH_MSG_DEBUG("   " << m_zBinsCustomLooping);
    ATH_MSG_DEBUG("   " << m_useVariableMiddleSPRange);
    if ( m_useVariableMiddleSPRange ) {
      ATH_MSG_DEBUG("   " << m_deltaRMiddleMinSPRange);
      ATH_MSG_DEBUG("   " << m_deltaRMiddleMaxSPRange);
    } else if ( not m_rRangeMiddleSP.empty() )
      ATH_MSG_DEBUG("   " << m_rRangeMiddleSP);
    ATH_MSG_DEBUG("   " << m_seedConfirmation);
    if ( m_seedConfirmation ) {
      ATH_MSG_DEBUG("   " << m_seedConfCentralZMin);
      ATH_MSG_DEBUG("   " << m_seedConfCentralZMax);
      ATH_MSG_DEBUG("   " << m_seedConfCentralRMax);
      ATH_MSG_DEBUG("   " << m_seedConfCentralNTopLargeR);
      ATH_MSG_DEBUG("   " << m_seedConfCentralNTopSmallR);
      ATH_MSG_DEBUG("   " << m_seedConfCentralMinBottomRadius);
      ATH_MSG_DEBUG("   " << m_seedConfCentralMaxZOrigin);
      ATH_MSG_DEBUG("   " << m_seedConfCentralMinImpact);
      ATH_MSG_DEBUG("   " << m_seedConfForwardZMin);
      ATH_MSG_DEBUG("   " << m_seedConfForwardZMax);
      ATH_MSG_DEBUG("   " << m_seedConfForwardRMax);
      ATH_MSG_DEBUG("   " << m_seedConfForwardNTopLargeR);
      ATH_MSG_DEBUG("   " << m_seedConfForwardNTopSmallR);
      ATH_MSG_DEBUG("   " << m_seedConfForwardMinBottomRadius);
      ATH_MSG_DEBUG("   " << m_seedConfForwardMaxZOrigin);
      ATH_MSG_DEBUG("   " << m_seedConfForwardMinImpact);
    }
    ATH_MSG_DEBUG("   " << m_useDetailedDoubleMeasurementInfo);
    ATH_MSG_DEBUG("   " << m_toleranceParam);
    ATH_MSG_DEBUG("   " << m_phiMin);
    ATH_MSG_DEBUG("   " << m_phiMax);
    ATH_MSG_DEBUG("   " << m_rMin);
    ATH_MSG_DEBUG("   " << m_zAlign);
    ATH_MSG_DEBUG("   " << m_rAlign);
    ATH_MSG_DEBUG("   " << m_sigmaError);

    ATH_MSG_DEBUG(" * Used by SeedFilterConfig:");
    ATH_MSG_DEBUG("   " << m_deltaRMin);
    ATH_MSG_DEBUG("   " << m_maxSeedsPerSpM);
    ATH_MSG_DEBUG("   " << m_useDeltaRorTopRadius);
    ATH_MSG_DEBUG("   " << m_curvatureSortingInFilter);
    ATH_MSG_DEBUG("   " << m_seedConfirmationInFilter);
    if (m_seedConfirmationInFilter) {
      ATH_MSG_DEBUG("   " << m_maxSeedsPerSpMConf);
      ATH_MSG_DEBUG("   " << m_maxQualitySeedsPerSpMConf);
      ATH_MSG_DEBUG("   " << m_seedConfCentralZMin);
      ATH_MSG_DEBUG("   " << m_seedConfCentralZMax);
      ATH_MSG_DEBUG("   " << m_seedConfCentralRMax);
      ATH_MSG_DEBUG("   " << m_seedConfCentralNTopLargeR);
      ATH_MSG_DEBUG("   " << m_seedConfCentralNTopSmallR);
      ATH_MSG_DEBUG("   " << m_seedConfCentralMinBottomRadius);
      ATH_MSG_DEBUG("   " << m_seedConfCentralMaxZOrigin);
      ATH_MSG_DEBUG("   " << m_seedConfCentralMinImpact);
      ATH_MSG_DEBUG("   " << m_seedConfForwardZMin);
      ATH_MSG_DEBUG("   " << m_seedConfForwardZMax);
      ATH_MSG_DEBUG("   " << m_seedConfForwardRMax);
      ATH_MSG_DEBUG("   " << m_seedConfForwardNTopLargeR);
      ATH_MSG_DEBUG("   " << m_seedConfForwardNTopSmallR);
      ATH_MSG_DEBUG("   " << m_seedConfForwardMinBottomRadius);
      ATH_MSG_DEBUG("   " << m_seedConfForwardMaxZOrigin);
      ATH_MSG_DEBUG("   " << m_seedConfForwardMinImpact);
    }
    ATH_MSG_DEBUG("   " << m_impactWeightFactor);
    ATH_MSG_DEBUG("   " << m_compatSeedWeight);
    ATH_MSG_DEBUG("   " << m_compatSeedLimit);
    ATH_MSG_DEBUG("   " << m_seedWeightIncrement);
    ATH_MSG_DEBUG("   " << m_numSeedIncrement);
    ATH_MSG_DEBUG("   " << m_deltaInvHelixDiameter);

 
    if (m_zBinEdges.size() - 1 !=
      m_zBinNeighborsTop.size() and
      not m_zBinNeighborsTop.empty()) {
      ATH_MSG_ERROR("Inconsistent config zBinNeighborsTop");
      return StatusCode::FAILURE;
    }

    if (m_zBinEdges.size() - 1 !=
      m_zBinNeighborsBottom.size() and
      not m_zBinNeighborsBottom.empty()) {
      ATH_MSG_ERROR("Inconsistent config zBinNeighborsBottom");
      return StatusCode::FAILURE;
    }

    if (m_zBinsCustomLooping.size() != 0) {
      // check if zBinsCustomLooping contains numbers from 1 to the total number
      // of bin in zBinEdges
      for (size_t i = 1; i != m_zBinEdges.size(); i++) {
        if (std::find(m_zBinsCustomLooping.begin(),
                      m_zBinsCustomLooping.end(),
                      i) == m_zBinsCustomLooping.end()) {
          ATH_MSG_ERROR("Inconsistent config zBinsCustomLooping does not contain the same bins as zBinEdges");
          return StatusCode::FAILURE;
        }
      }
    }

    return StatusCode::SUCCESS;
  }

  StatusCode
  SeedingTool::createSeeds(const EventContext& /*ctx*/,
  		           const std::vector<const ActsTrk::SpacePoint*>& spContainer,
			   const Acts::Vector3& beamSpotPos,
			   const Acts::Vector3& bField,
			   ActsTrk::SeedContainer& seedContainer ) const
  {					 
    // Create Seeds
    //TODO POSSIBLE OPTIMISATION come back here: see MR !52399 ( i.e. use static thread_local)
    std::vector< seed_type > groupSeeds;
    ATH_CHECK(createSeeds(spContainer.begin(),
			  spContainer.end(),
			  beamSpotPos,
			  bField,
			  groupSeeds));
    
    // Store seeds
    seedContainer.reserve(groupSeeds.size());
    for( const auto& seed: groupSeeds) {
      std::unique_ptr< seed_type > to_add = 
	std::make_unique< seed_type >(seed);
      seedContainer.push_back(std::move(to_add));  
    }
    
    return StatusCode::SUCCESS;
  }

  template< typename external_iterator_t >
  StatusCode
  SeedingTool::createSeeds(external_iterator_t spBegin,
			   external_iterator_t spEnd,
			   const Acts::Vector3& beamSpotPos,
			   const Acts::Vector3& bField,
			   std::vector< seed_type >& seeds) const 
  {
    seeds.clear();
    if (spBegin == spEnd)
      return StatusCode::SUCCESS;
    
    auto [gridCfg, finderCfg] = prepareConfiguration(Acts::Vector2(beamSpotPos[Amg::x], beamSpotPos[Amg::y]),
    	 	   	      			     bField); 
        
    auto extractCovariance = [&beamSpotPos] (const value_type& sp, 
					     float, float, float) -> std::pair<Acts::Vector3, Acts::Vector2> {
      /// Convert coordinates w.r.t. beam spot
      Acts::Vector3 position(sp.x() - beamSpotPos[Amg::x], 
			     sp.y() - beamSpotPos[Amg::y], 
			     sp.z() - beamSpotPos[Amg::z]);
      Acts::Vector2 covariance(sp.varianceR(), sp.varianceZ());
      return std::make_pair(position, covariance);
    };
    
    
    Acts::Extent rRangeSPExtent;
    
    std::shared_ptr< Acts::BinFinder< value_type > > bottomBinFinder =
      std::make_shared< Acts::BinFinder< value_type > >(m_zBinNeighborsBottom, m_numPhiNeighbors);
    std::shared_ptr< Acts::BinFinder< value_type > > topBinFinder =
      std::make_shared< Acts::BinFinder< value_type > >(m_zBinNeighborsTop, m_numPhiNeighbors);
    
    std::unique_ptr< Acts::SpacePointGrid< value_type > > grid =
      Acts::SpacePointGridCreator::createGrid< value_type >(gridCfg);
    Acts::BinnedSPGroup< value_type > spacePointsGrouping(spBegin, spEnd, extractCovariance,								     
      bottomBinFinder, topBinFinder, std::move(grid), rRangeSPExtent, finderCfg);
    
    Acts::SeedFinder< value_type > finder(finderCfg);

    // variable middle SP radial region of interest
    const Acts::Range1D<float> rMiddleSPRange(std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
					      finderCfg.deltaRMiddleMinSPRange,
					      std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2 -
					      finderCfg.deltaRMiddleMaxSPRange);

    //TODO POSSIBLE OPTIMISATION come back here: see MR !52399 ( i.e. use static thread_local)
    typename decltype(finder)::SeedingState state;

    auto group = spacePointsGrouping.begin();
    auto groupEnd = spacePointsGrouping.end();
    for (; group != groupEnd; ++group) {
      finder.createSeedsForGroup(state, std::back_inserter(seeds), group.bottom(),
				 group.middle(), group.top(), rMiddleSPRange);
    }

    return StatusCode::SUCCESS;
  }


  const std::pair< Acts::SpacePointGridConfig, 
		   Acts::SeedFinderConfig< typename SeedingTool::value_type > >
  SeedingTool::prepareConfiguration(const Acts::Vector2& beamPos,
				    const Acts::Vector3& bField) const {
    
    //TODO POSSIBLE OPTIMISATION
    // do not create for each call SpacePointGridConfig, SeedFinderConfig and SeedFilterConfig.
    // They only two quantities that depend on the event context are
    // finderCfg.bFieldInZ and finderCfg.beamPos.
    // The rest of the configuration stays the same.
    
    // Configuration for Acts::SpacePointGrid
    // These values will not be changed during execution.
    // For the grid formation, a constant value of the Z component
    // is used to reduce computational overhead.
    Acts::SpacePointGridConfig gridCfg;
    gridCfg.minPt = m_minPt;
    gridCfg.cotThetaMax = m_cotThetaMax;
    gridCfg.impactMax = m_impactMax;
    gridCfg.zMin = m_zMin;
    gridCfg.zMax = m_zMax;
    gridCfg.phiMin = m_gridPhiMin;
    gridCfg.phiMax = m_gridPhiMax;
    gridCfg.zBinEdges = m_zBinEdges;
    gridCfg.deltaRMax = m_deltaRMax;
    gridCfg.rMax = m_gridRMax;
    gridCfg.bFieldInZ = m_bFieldInZ;
    gridCfg.phiBinDeflectionCoverage = m_phiBinDeflectionCoverage;

    // Configuration for Acts::SeedFinder
    // These values will not be changed during execution
    // B Field and Beam Spot position will be updated for each event (finderCfg.bFieldInZ and finderCfg.beamPos)
    Acts::SeedFinderConfig< value_type > finderCfg;
    finderCfg.bFieldInZ = bField[2];
    finderCfg.beamPos = beamPos;
    finderCfg.minPt = m_minPt;
    finderCfg.cotThetaMax = m_cotThetaMax;
    finderCfg.impactMax = m_impactMax;
    finderCfg.zMin = m_zMin;
    finderCfg.zMax = m_zMax;
    finderCfg.zBinEdges = m_zBinEdges;
    finderCfg.rMax = m_rMax;
    finderCfg.binSizeR = m_binSizeR;
    finderCfg.forceRadialSorting = m_forceRadialSorting;
    finderCfg.deltaRMin = m_deltaRMin;
    finderCfg.deltaRMax = m_deltaRMax;
    finderCfg.deltaRMinTopSP = m_deltaRMinTopSP;
    finderCfg.deltaRMaxTopSP = m_deltaRMaxTopSP;
    finderCfg.deltaRMinBottomSP = m_deltaRMinBottomSP;
    finderCfg.deltaRMaxBottomSP = m_deltaRMaxBottomSP;
    finderCfg.deltaZMax = m_deltaZMax;
    finderCfg.collisionRegionMin = m_collisionRegionMin;
    finderCfg.collisionRegionMax = m_collisionRegionMax;
    finderCfg.sigmaScattering = m_sigmaScattering;
    finderCfg.maxPtScattering = m_maxPtScattering;
    finderCfg.radLengthPerSeed = m_radLengthPerSeed;
    finderCfg.maxSeedsPerSpM = m_maxSeedsPerSpM;
    finderCfg.interactionPointCut = m_interactionPointCut;
    finderCfg.arithmeticAverageCotTheta = m_arithmeticAverageCotTheta;
    finderCfg.skipPreviousTopSP = m_skipPreviousTopSP;
    finderCfg.zBinsCustomLooping = m_zBinsCustomLooping;
    finderCfg.useVariableMiddleSPRange = m_useVariableMiddleSPRange;
    finderCfg.deltaRMiddleMinSPRange = m_deltaRMiddleMinSPRange;
    finderCfg.deltaRMiddleMaxSPRange = m_deltaRMiddleMaxSPRange;
    finderCfg.seedConfirmation = m_seedConfirmation;
    finderCfg.centralSeedConfirmationRange.zMinSeedConf = m_seedConfCentralZMin;
    finderCfg.centralSeedConfirmationRange.zMaxSeedConf = m_seedConfCentralZMax;
    finderCfg.centralSeedConfirmationRange.rMaxSeedConf = m_seedConfCentralRMax;
    finderCfg.centralSeedConfirmationRange.nTopForLargeR = m_seedConfCentralNTopLargeR;
    finderCfg.centralSeedConfirmationRange.nTopForSmallR = m_seedConfCentralNTopSmallR;
    finderCfg.centralSeedConfirmationRange.seedConfMinBottomRadius = m_seedConfCentralMinBottomRadius;
    finderCfg.centralSeedConfirmationRange.seedConfMaxZOrigin = m_seedConfCentralMaxZOrigin;
    finderCfg.centralSeedConfirmationRange.minImpactSeedConf = m_seedConfCentralMinImpact;
    finderCfg.forwardSeedConfirmationRange.zMinSeedConf = m_seedConfForwardZMin;
    finderCfg.forwardSeedConfirmationRange.zMaxSeedConf = m_seedConfForwardZMax;
    finderCfg.forwardSeedConfirmationRange.rMaxSeedConf = m_seedConfForwardRMax;
    finderCfg.forwardSeedConfirmationRange.nTopForLargeR = m_seedConfForwardNTopLargeR;
    finderCfg.forwardSeedConfirmationRange.nTopForSmallR = m_seedConfForwardNTopSmallR;
    finderCfg.forwardSeedConfirmationRange.seedConfMinBottomRadius = m_seedConfForwardMinBottomRadius;
    finderCfg.forwardSeedConfirmationRange.seedConfMaxZOrigin = m_seedConfForwardMaxZOrigin;
    finderCfg.forwardSeedConfirmationRange.minImpactSeedConf = m_seedConfForwardMinImpact;
    finderCfg.useDetailedDoubleMeasurementInfo = m_useDetailedDoubleMeasurementInfo;
    finderCfg.toleranceParam = m_toleranceParam;
    finderCfg.phiMin = m_phiMin;
    finderCfg.phiMax = m_phiMax;
    finderCfg.rMin = m_rMin;
    finderCfg.zAlign = m_zAlign;
    finderCfg.rAlign = m_rAlign;
    finderCfg.sigmaError = m_sigmaError;

    if (m_useDetailedDoubleMeasurementInfo) {
      finderCfg.getTopHalfStripLength.connect(
        [](const void*, const value_type& sp) -> float {
          return sp.topHalfStripLength();
        });
      finderCfg.getBottomHalfStripLength.connect(
        [](const void*, const value_type& sp) -> float {
          return sp.bottomHalfStripLength();
        });
      finderCfg.getTopStripDirection.connect(
        [](const void*, const value_type& sp) -> Acts::Vector3 {
          return sp.topStripDirection();
        });
      finderCfg.getBottomStripDirection.connect(
        [](const void*, const value_type& sp) -> Acts::Vector3 {
          return sp.bottomStripDirection();
        });
      finderCfg.getStripCenterDistance.connect(
          [](const void*, const value_type& sp) -> Acts::Vector3 {
            return sp.stripCenterDistance();
          });
      finderCfg.getTopStripCenterPosition.connect(
          [](const void*, const value_type& sp) -> Acts::Vector3 {
            return sp.topStripCenter();
          });
    }

    // Configuration for Acts::SeedFilter
    Acts::SeedFilterConfig filterCfg;
    filterCfg.deltaRMin = m_deltaRMin;
    filterCfg.maxSeedsPerSpM = m_maxSeedsPerSpM;
    filterCfg.useDeltaRorTopRadius = m_useDeltaRorTopRadius;
    filterCfg.curvatureSortingInFilter = m_curvatureSortingInFilter;
    filterCfg.seedConfirmation = m_seedConfirmationInFilter;
    filterCfg.maxSeedsPerSpMConf = m_maxSeedsPerSpMConf;
    filterCfg.maxQualitySeedsPerSpMConf = m_maxQualitySeedsPerSpMConf;
    filterCfg.centralSeedConfirmationRange = finderCfg.centralSeedConfirmationRange;
    filterCfg.forwardSeedConfirmationRange = finderCfg.forwardSeedConfirmationRange;
    filterCfg.impactWeightFactor = m_impactWeightFactor;
    filterCfg.compatSeedWeight = m_compatSeedWeight;
    filterCfg.compatSeedLimit = m_compatSeedLimit;
    filterCfg.seedWeightIncrement = m_seedWeightIncrement;
    filterCfg.numSeedIncrement = m_numSeedIncrement;
    filterCfg.deltaInvHelixDiameter = m_deltaInvHelixDiameter;
    finderCfg.seedFilter = std::make_unique<Acts::SeedFilter< value_type > >(filterCfg);

    return std::make_pair(gridCfg, finderCfg);
  }

} // namespace ActsTrk
