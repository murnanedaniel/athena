/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#include "MuonHoughPatternTools/MuonLayerHoughTool.h"
#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "MuonReadoutGeometry/sTgcReadoutElement.h"
#include "MuonReadoutGeometry/MuonChannelDesign.h"
#include "MuonReadoutGeometry/MuonPadDesign.h"
#include "MuonPattern/MuonPatternCombination.h"
#include "MuonPattern/MuonPatternChamberIntersect.h"
#include "TFile.h"
#include "TTree.h"
#include "AtlasHepMC/GenEvent.h"
#include "GaudiKernel/ConcurrencyFlags.h"
#include "CxxUtils/sincos.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODMuon/MuonSegmentContainer.h"
namespace Muon {

  MuonLayerHoughTool::MuonLayerHoughTool(const std::string& type, const std::string& name, const IInterface* parent):
    AthAlgTool(type,name,parent),
    m_ntechnologies(UINT_MAX), // gets set inside initialize()
    m_techToTruthNameIdx()
  {
    declareInterface<MuonLayerHoughTool>(this);
    declareInterface<IMuonHoughPatternFinderTool>(this);
  }

  StatusCode MuonLayerHoughTool::initialize() {

    ATH_CHECK( m_idHelperSvc.retrieve() );
    m_ntechnologies = m_idHelperSvc->mdtIdHelper().technologyNameIndexMax()+1;
    ATH_CHECK( m_printer.retrieve() );
    if( m_doTruth && !m_truthSummaryTool.empty() ){
      ATH_CHECK( m_truthSummaryTool.retrieve() );
    } else {
      m_truthSummaryTool.disable();
    }
    const MuonGM::MuonDetectorManager* muDetMgr=nullptr;
    ATH_CHECK( detStore()->retrieve( muDetMgr ) );

    if( m_doNtuple ){
      if (Gaudi::Concurrency::ConcurrencyFlags::concurrent() && Gaudi::Concurrency::ConcurrencyFlags::numThreads()>1) {
        // Disabled for >1 threads due to thread-safety concerns, but we want to keep it as a debug tool
        ATH_MSG_DEBUG("HitNtuple disabled because of concurrency");
	m_file = 0;
	m_tree = 0;
	m_ntuple = 0;
      } else {
        TDirectory* cdir = gDirectory;
        m_file = new TFile("HitNtuple.root","RECREATE");
        m_tree = new TTree("data","data");
        m_ntuple = new MuonHough::HitNtuple();
        m_ntuple->initForWrite(*m_tree);
        gDirectory = cdir;
      }
    }else{
      m_file = 0;
      m_tree = 0;
      m_ntuple = 0;
    }
    
    initializeSectorMapping(muDetMgr);

    // if m_truthNames is empty, fill it if running on truth
    if( m_truthNames.empty() && m_doTruth ){
      std::string postfix = "_TruthMap";
      std::string allNames("");
      unsigned int i=0;
      for( unsigned int tech=0; tech<m_ntechnologies;++tech ){
        // check if technology is part of layout
        if (tech==MuonStationIndex::CSC && !m_idHelperSvc->hasCSC()) continue;
        else if (tech==MuonStationIndex::STGC && !m_idHelperSvc->hasSTgc()) continue;
        else if (tech==MuonStationIndex::MM && !m_idHelperSvc->hasMM()) continue;
        std::string thisname = std::string(m_idHelperSvc->mdtIdHelper().technologyString(tech)) + postfix;
        m_truthNames.emplace_back( thisname );
        // since we need to access the elements of m_truthNames later on, we need to remember 
        // which technology is saved at which index of the vector
        m_techToTruthNameIdx.insert(std::make_pair(tech, i));
        ++i;
        allNames += " ";
        allNames += thisname;
      }
      ATH_MSG_DEBUG("TruthMaps " << allNames );
    }
    if(!m_doTruth){ m_truthNames.clear(); } //Nullify if not using collections

    ATH_CHECK( m_truthNames.initialize() );
    if(m_doNtuple && m_doTruth){
      ATH_CHECK(m_MuonTruthParticlesKey.initialize());
      ATH_CHECK(m_MuonTruthSegmentsKey.initialize());
    }
    else{
      m_MuonTruthParticlesKey="";
      m_MuonTruthSegmentsKey="";
    }

    // initialize cuts, if only one cut, use make_pair to avoid compiler issues, format is (position, cut)
    m_selectors.resize(MuonStationIndex::ChIndexMax);
    m_selectors[MuonStationIndex::BIS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,5.9)}); // old values: 6.9; optimized: 7.9
    m_selectors[MuonStationIndex::BIL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,5.9)}); // old values: 6.9; optimized: 7.9
    m_selectors[MuonStationIndex::BMS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 7.9; optimized: 7.9
    m_selectors[MuonStationIndex::BML] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 7.9; optimized: 7.9 
    m_selectors[MuonStationIndex::BOS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 4.9; optimized: 5.9
    m_selectors[MuonStationIndex::BOL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 4.9; optimized: 5.9
    m_selectors[MuonStationIndex::BEE] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,3.9)}); // old values: 5.9; optimized: 5.9
    m_selectors[MuonStationIndex::EIS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,5.9)}); // old values: 6.9; optimized: 7.9
    m_selectors[MuonStationIndex::EIL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,5.9)}); // old values: 6.9; optimized: 7.9
    m_selectors[MuonStationIndex::EMS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 7.9; optimized: 5.9
    m_selectors[MuonStationIndex::EML] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 7.9; optimized: 5.9
    m_selectors[MuonStationIndex::EOS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 4.9; optimized: 5.9
    m_selectors[MuonStationIndex::EOL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 4.9; optimized: 5.9
    m_selectors[MuonStationIndex::EES] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 4.9; optimized: 5.9
    m_selectors[MuonStationIndex::EEL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,4.9)}); // old values: 4.9; optimized: 5.9

    m_selectorsLoose.resize(MuonStationIndex::ChIndexMax);
    m_selectorsLoose[MuonStationIndex::BIS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,1.9)}); // old values: 2.9; optimized: 3.9
    m_selectorsLoose[MuonStationIndex::BIL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 2.9; optimized: 3.9
    m_selectorsLoose[MuonStationIndex::BMS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,1.9)}); // old values: 4.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::BML] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 4.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::BOS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,1.9)}); // old values: 2.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::BOL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 2.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::BEE] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,1.9)}); // old values: 3.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::EIS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,1.9)}); // old values: 4.9; optimized: 3.9
    m_selectorsLoose[MuonStationIndex::EIL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,1.9)}); // old values: 4.9; optimized: 3.9
    m_selectorsLoose[MuonStationIndex::EMS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 5.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::EML] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 5.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::EOS] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 2.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::EOL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 2.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::EES] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 2.9; optimized: 2.9
    m_selectorsLoose[MuonStationIndex::EEL] = MuonHough::MuonLayerHoughSelector({std::make_pair(0,2.9)}); // old values: 2.9; optimized: 2.9

    return StatusCode::SUCCESS;
  }

  StatusCode MuonLayerHoughTool::finalize() {
    // ??? finalizeNTuple -> may be moved somewhere
    if( m_doNtuple && m_ntuple){
      TDirectory* cdir = gDirectory;
      m_file->cd();
      m_tree->Write();
      m_file->Write();
      m_file->Close();
      delete m_ntuple;
      gDirectory = cdir;
    }
    return StatusCode::SUCCESS;
  }

  void MuonLayerHoughTool::getTruth() const {
    if(m_ntuple){
      SG::ReadHandle<xAOD::TruthParticleContainer> truthMuons(m_MuonTruthParticlesKey);
      if (truthMuons.isValid()) {
        ATH_MSG_DEBUG("Retrieved truth muons " << truthMuons->size());
        int nmuons = 0;
        for (const auto truthMu: *truthMuons){
          m_ntuple->tpdgId[nmuons] = truthMu->pdgId();
          m_ntuple->tbarcode[nmuons] = truthMu->barcode();
          m_ntuple->tmuonIndex[nmuons] = nmuons;
          m_ntuple->pt[nmuons] = truthMu->pt();
          m_ntuple->eta[nmuons] = truthMu->eta();
          m_ntuple->phi[nmuons] = truthMu->phi();
          m_ntuple->nmdts[nmuons] = 0;
          m_ntuple->nrpcs[nmuons] = 0;
          m_ntuple->ntgcs[nmuons] = 0;
          m_ntuple->ncscs[nmuons] = 0;
          m_ntuple->ntmdts[nmuons] = 0;
          m_ntuple->ntrpcs[nmuons] = 0;
          m_ntuple->nttgcs[nmuons] = 0;
          m_ntuple->ntcscs[nmuons] = 0;
          ++nmuons;
        }
        m_ntuple->nmuons = nmuons;
        SG::ReadHandle<xAOD::MuonSegmentContainer> truthSegments(m_MuonTruthSegmentsKey);
        if (truthSegments.isValid()) {
          ATH_MSG_DEBUG("Retrieved truth Segments " << truthSegments->size());
          int nsegs = 0;
          for (const auto truthSeg: *truthSegments){
            m_ntuple->sbarcode[nsegs] = 0;
            m_ntuple->sposx[nsegs] = truthSeg->x();
            m_ntuple->sposy[nsegs] = truthSeg->y();
            m_ntuple->sposz[nsegs] = truthSeg->z();
            m_ntuple->sdirx[nsegs] = truthSeg->px();
            m_ntuple->sdiry[nsegs] = truthSeg->py();
            m_ntuple->sdirz[nsegs] = truthSeg->pz();
            m_ntuple->snPrecHits[nsegs] = truthSeg->nPrecisionHits();
            m_ntuple->snTrigHits[nsegs] = truthSeg->nPhiLayers() + truthSeg->nTrigEtaLayers();
            m_ntuple->sSector[nsegs] = truthSeg->sector();
            m_ntuple->sChIndex[nsegs] = truthSeg->chamberIndex();
            ++nsegs;
          }
          m_ntuple->nsegs = nsegs;
        }
      }
    }

  }
  
  void MuonLayerHoughTool::reset() const {
    if( m_ntuple )  m_ntuple->reset();
  }


  auto MuonLayerHoughTool::find( 
      const std::vector<const MdtPrepDataCollection*>& mdtCols,  
      const std::vector<const CscPrepDataCollection*>& ,  
      const std::vector<const TgcPrepDataCollection*>& tgcCols,  
      const std::vector<const RpcPrepDataCollection*>& rpcCols,  
      const MuonSegmentCombinationCollection* ) const 
      -> std::pair<std::unique_ptr<MuonPatternCombinationCollection>, std::unique_ptr<HoughDataPerSectorVec>> {
    reset();
    State state;
    ATH_MSG_DEBUG("MuonLayerHoughTool::find");
    if( m_doTruth ) getTruth();


    // create structure to hold data per sector and set the sector indices
    state.houghDataPerSectorVec->vec.resize(16);
    for( unsigned int i=0;i<state.houghDataPerSectorVec->vec.size();++i ) state.houghDataPerSectorVec->vec[i].sector=i+1;

    // return DetectorRegionIndex and sectorLayerHash
    auto getHashes = [this]( const Identifier& id ){
      MuonStationIndex::DetectorRegionIndex regionIndex = m_idHelperSvc->regionIndex(id);
      MuonStationIndex::LayerIndex layerIndex = m_idHelperSvc->layerIndex(id);
      unsigned int sectorLayerHash = MuonStationIndex::sectorLayerHash(regionIndex,layerIndex);
      return std::make_pair(regionIndex,sectorLayerHash);
    };

    for( auto col : mdtCols ){
      if( !col ) continue;
      Identifier id = col->identify();
      int sector = m_idHelperSvc->sector(id);
      auto hashes = getHashes(id);
      fill(state.truthHits,*col,state.houghDataPerSectorVec->vec[sector-1].hitVec[hashes.second]);
    }

    for( auto col : rpcCols ){
      if( !col ) continue;
      Identifier id = col->identify();
      int sector = m_idHelperSvc->sector(id);
      auto hashes = getHashes(id);
      fill(state.truthHits,*col,state.houghDataPerSectorVec->vec[sector-1].hitVec[hashes.second],state.houghDataPerSectorVec->vec[sector-1].phiHitVec[hashes.first]);
    }

    auto hashInSector = [this]( IdentifierHash hash, int sector, unsigned int sectorLayerHash ) {
      const std::vector<IdentifierHash>& hashes = m_collectionsPerSector[sector-1].technologyRegionHashVecs[MuonStationIndex::TGC][sectorLayerHash];
      return std::binary_search(hashes.begin(),hashes.end(),hash);
    };

    for( auto col : tgcCols ){
      if( !col ) continue;
      Identifier id = col->identify();
      int sector = m_idHelperSvc->sector(id);
      auto hashes = getHashes(id);
      // fill current sector
      fill(state.truthHits, state.houghDataPerSectorVec->tgcClusteringObjs, *col,state.houghDataPerSectorVec->vec[sector-1].hitVec[hashes.second],
          state.houghDataPerSectorVec->vec[sector-1].phiHitVec[hashes.first],sector);

      // fill neighbours if in overlap
      int neighbourSectorDown = sector == 1 ? 16 : sector-1;
      if( hashInSector(col->identifyHash(),neighbourSectorDown,hashes.second) ) 
        fill(state.truthHits, state.houghDataPerSectorVec->tgcClusteringObjs, *col,state.houghDataPerSectorVec->vec[neighbourSectorDown-1].hitVec[hashes.second],
             state.houghDataPerSectorVec->vec[neighbourSectorDown-1].phiHitVec[hashes.first],neighbourSectorDown);

      int neighbourSectorUp   = sector == 16 ? 1 : sector+1;
      if( hashInSector(col->identifyHash(),neighbourSectorUp,hashes.second) ) 
        fill(state.truthHits, state.houghDataPerSectorVec->tgcClusteringObjs, *col,state.houghDataPerSectorVec->vec[neighbourSectorUp-1].hitVec[hashes.second],
             state.houghDataPerSectorVec->vec[neighbourSectorUp-1].phiHitVec[hashes.first],neighbourSectorUp);
      
    }
    
    return analyse(state);
  }

  auto MuonLayerHoughTool::analyse(
      const MdtPrepDataContainer*  mdtCont,
      const CscPrepDataContainer*  cscCont,
      const TgcPrepDataContainer*  tgcCont,
      const RpcPrepDataContainer*  rpcCont,
      const sTgcPrepDataContainer* stgcCont,  
      const MMPrepDataContainer*   mmCont ) const 
      -> std::pair<std::unique_ptr<MuonPatternCombinationCollection>, std::unique_ptr<HoughDataPerSectorVec>> {
    reset();
    State state; 
    ATH_MSG_DEBUG("MuonLayerHoughTool::analyse");
    if( m_doTruth ) getTruth();

    state.houghDataPerSectorVec->vec.resize(16);

    // loops over all sectors, contains hashes for technology and region and chamber (?)
    CollectionsPerSectorCit sit = m_collectionsPerSector.begin();
    CollectionsPerSectorCit sit_end = m_collectionsPerSector.end();
    for( ;sit!=sit_end;++sit){
            
      ATH_MSG_DEBUG("analyse: Filling hits sector " << sit->sector);

      HoughDataPerSector& houghData = state.houghDataPerSectorVec->vec[sit->sector-1]; 
      houghData.sector = sit->sector;

      // fill hits for this sector -> hitsVec and PhiHitsVec are known now
      fillHitsPerSector( state.truthHits, state.houghDataPerSectorVec->tgcClusteringObjs, *sit, mdtCont,cscCont,tgcCont,rpcCont,stgcCont,mmCont,houghData);

    }
    return analyse(state);
  }
  
  auto MuonLayerHoughTool::analyse(State& state) const 
      -> std::pair<std::unique_ptr<MuonPatternCombinationCollection>, std::unique_ptr<HoughDataPerSectorVec>> {    

    auto patternCombis = std::make_unique<MuonPatternCombinationCollection>();
    
    // loop over data and fill the hough transform
    for( auto& houghData : state.houghDataPerSectorVec->vec ){

      ATH_MSG_DEBUG("analyse: Filling Hough sector " << houghData.sector);

      // loop over all possible station layers in the sector and run the eta transform
      for( unsigned int layerHash=0;layerHash<MuonStationIndex::sectorLayerHashMax();++layerHash ){
        
        // get hits for layer, skip empty layers
        HitVec& hits = houghData.hitVec[layerHash];
        if( hits.empty() ) continue;

        // decompose hash, calculate indices etc
        auto regionLayer = MuonStationIndex::decomposeSectorLayerHash( layerHash );
        MuonStationIndex::DetectorRegionIndex region = regionLayer.first;
        MuonStationIndex::LayerIndex          layer  = regionLayer.second;
        MuonStationIndex::StIndex             index = MuonStationIndex::toStationIndex(region,layer);

        // get Hough transform
        MuonHough::MuonLayerHough& hough = state.houghDataPerSectorVec->detectorHoughTransforms.hough( houghData.sector,  region, layer );
       
        ATH_MSG_DEBUG("analyse: Filling Summary: loc s"   << houghData.sector 
                      << " "           << MuonStationIndex::regionName(region) 
                      << " "           << MuonStationIndex::layerName(layer) 
                      << " -> stIndex: " << MuonStationIndex::stName(index) 
                      << " etaHits: "  << hits.size() );

        // look for maxima using hough in eta per layer
        if( !findMaxima(state.truthHits, state.foundTruthHits, state.seedMaxima, hough,hits,
          houghData.maxVec[layerHash]) || houghData.maxVec[layerHash].empty() ) continue;

        ++houghData.nlayersWithMaxima[region];
        houghData.nmaxHitsInRegion[region] += houghData.maxVec[layerHash].front()->max;

        ATH_MSG_DEBUG("analyse: Eta maxima Summary: loc s"   << houghData.sector 
                      << " "           << MuonStationIndex::regionName(region) 
                      << " "           << MuonStationIndex::layerName(layer) 
                      << " -> stIndex: " << MuonStationIndex::stName(index) 
                      << " hash: "     << layerHash 
                      << " nMaxima: "  << houghData.maxVec[layerHash].size() );
      }  // loop over layerHash -> maxima per layer in eta are known now
    } // loop over sectors

    if( m_useSeeds ){
      std::vector<Road> roads;
      buildRoads(state.seedMaxima, state.houghDataPerSectorVec->detectorHoughTransforms, state.truthHits, state.foundTruthHits, 
        state.houghDataPerSectorVec, roads);
      
      // create association map
      ATH_MSG_DEBUG("analyse: Building pattern combinations using roads " << roads.size() );
      for( auto& road : roads ){
        std::map< MuonHough::MuonPhiLayerHough::Maximum*, MuonLayerHoughTool::MaximumVec > phiEtaAssMap;
        MuonLayerHoughTool::RegionMaximumVec                                               unassociatedEtaMaxima;

        int                                   sector  = road.seed->hough->m_descriptor.sector;
        MuonStationIndex::ChIndex             chIndex = road.seed->hough->m_descriptor.chIndex;
        MuonStationIndex::LayerIndex          layer   = Muon::MuonStationIndex::toLayerIndex(chIndex);
        MuonStationIndex::DetectorRegionIndex region  = road.seed->hough->m_descriptor.region;
        ATH_MSG_DEBUG("analyse: Seeding new road: eta maxima " << road.maxima.size() << " phi " << road.phiMaxima.size()
                      << " seed : sector " << sector << " " << Muon::MuonStationIndex::regionName(region) 
                      << " " << Muon::MuonStationIndex::layerName(layer)
                      << " maximum " << road.seed->max << " position " << road.seed->pos << " angle " << road.seed->theta );

        if( road.phiMaxima.empty() ) unassociatedEtaMaxima.push_back(road.maxima);
        else{
          for( auto& max : road.mergedPhiMaxima ) {
            phiEtaAssMap[&max] = road.maxima;
          }
        }
        createPatternCombinations(state.truthHits,state.outputTruthHits,phiEtaAssMap,*patternCombis);
        createPatternCombinations(unassociatedEtaMaxima,*patternCombis);
      }

    }
    else{
      // now that the full hough transform is filled, order sectors by maxima
      std::vector<HoughDataPerSector*> sectorData(state.houghDataPerSectorVec->vec.size());
      for( unsigned int i=0;i<state.houghDataPerSectorVec->vec.size();++i) sectorData[i] = &state.houghDataPerSectorVec->vec[i];
      std::stable_sort(sectorData.begin(),sectorData.end(),SortHoughDataPerSector());

      std::vector<HoughDataPerSector*>::iterator spit = sectorData.begin();
      std::vector<HoughDataPerSector*>::iterator spit_end = sectorData.end();
      for( ;spit!=spit_end;++spit ){

        // get data for this sector
        HoughDataPerSector& houghData = **spit;

        // loop over regions
        for( int reg = 0;reg<MuonStationIndex::DetectorRegionIndexMax;++reg ){   

          MuonStationIndex::DetectorRegionIndex region = static_cast<MuonStationIndex::DetectorRegionIndex>(reg);
  
          // only run analysis on sectors with maxima
          if( houghData.nlayersWithMaxima[region] == 0 ) continue;
          ATH_MSG_DEBUG("Analyzing sector " << (*spit)->sector << " " << MuonStationIndex::regionName(region) << " nmax " << (*spit)->maxEtaHits() 
                        << " layers with eta maxima " << houghData.nlayersWithMaxima[region] << " hits " << houghData.nmaxHitsInRegion[region]
                        << " layers with phi maxima " << houghData.nphilayersWithMaxima[region] << " hits " << houghData.nphimaxHitsInRegion[region] );
  
          // look for maxima in the overlap regions of sectors
          associateMaximaInNeighbouringSectors(houghData,state.houghDataPerSectorVec->vec);

          // layers in this region
          int nlayers = MuonStationIndex::LayerIndexMax;

          // first link phi maxima with eta maxima
          RegionMaximumVec unassociatedEtaMaxima(nlayers);
          std::map< MuonHough::MuonPhiLayerHough::Maximum*, MaximumVec > phiEtaAssociations;
          associateMaximaToPhiMaxima( region, houghData,  phiEtaAssociations, unassociatedEtaMaxima );

          // create pattern combinations for combined patterns 
          createPatternCombinations(state.truthHits, state.outputTruthHits, phiEtaAssociations, *patternCombis );
  
          // create pattern combinations for unassociated patterns 
          createPatternCombinations(unassociatedEtaMaxima,*patternCombis);
        }
      }
    }
    
    if( m_ntuple ) {
      fillNtuple(*(state.houghDataPerSectorVec.get()));
      m_tree->Fill();
    }

    ATH_MSG_DEBUG("Found " << patternCombis->size() << " pattern combinations " << std::endl << m_printer->print( *patternCombis ));

    if( m_doTruth && msgLvl(MSG::DEBUG) ){
      ATH_MSG_DEBUG("Hough performance ");
      printTruthSummary(state.truthHits,state.foundTruthHits);
      ATH_MSG_DEBUG("Association performance ");
      printTruthSummary(state.foundTruthHits,state.outputTruthHits);
    }

    return {std::move(patternCombis), std::move(state.houghDataPerSectorVec)};
  }

  void MuonLayerHoughTool::buildRoads(MaximumVec& seedMaxima, MuonHough::MuonDetectorHough& detectorHoughTransforms,
      std::set<Identifier>& truthHits, std::set<Identifier>& foundTruthHits, 
      std::unique_ptr<HoughDataPerSectorVec>& houghDataPerSectorVec, std::vector<MuonLayerHoughTool::Road>& roads ) const {
    // sort maxima according to hits
    std::stable_sort( seedMaxima.begin(),seedMaxima.end(), [](const MuonHough::MuonLayerHough::Maximum* m1,
                                                                  const MuonHough::MuonLayerHough::Maximum* m2 ){ return m1->max > m2->max; } );
    // loop over seed maxima (which are maxima) that pass certain thresholds detailed in cut_values
    std::set<const MuonHough::MuonLayerHough::Maximum*> associatedMaxima;
    for( auto maxit = seedMaxima.begin();maxit!=seedMaxima.end(); ++maxit ){
      // if this maximum is already in the set of associated maxima, do not do anything
      if( associatedMaxima.count(*maxit) ) continue;

      // maximum becomes our new seed
      MuonHough::MuonLayerHough::Maximum& seed = **maxit;
      
      // decomposing the locality information for the seed
      int                                   sector  = seed.hough->m_descriptor.sector;
      MuonStationIndex::ChIndex             chIndex = seed.hough->m_descriptor.chIndex;
      MuonStationIndex::LayerIndex          layer   = Muon::MuonStationIndex::toLayerIndex(chIndex);
      MuonStationIndex::DetectorRegionIndex region  = seed.hough->m_descriptor.region;
      
      // creating new road with said seed
      Road road(seed);
      
      ATH_MSG_DEBUG(" New seed: sector " << seed.hough->m_descriptor.sector << " " << Muon::MuonStationIndex::regionName(region) 
                    << " " << Muon::MuonStationIndex::layerName(layer)
                    << " maximum " << seed.max << " position " << seed.pos << " angle " << seed.theta << " ptr " << &seed );

      bool isNSW=m_idHelperSvc->issTgc(seed.hits[0]->prd->identify()) || m_idHelperSvc->isMM(seed.hits[0]->prd->identify());
      // extend seed within the current sector
      // sector indices have an offset of -1 because the numbering of the sectors are from 1 to 16 but the indices in the vertices are of course 0 to 15
      extendSeed( detectorHoughTransforms, truthHits, foundTruthHits, road, houghDataPerSectorVec->vec[sector-1] );

      // look for maxima in the overlap regions of sectors
      int sectorN            = sector-1;
      if(sectorN<1)  sectorN = 16;
      int sectorP            = sector+1;
      if(sectorP>16) sectorP = 1;

      // associate the road with phi maxima
      associatePhiMaxima( road, houghDataPerSectorVec->vec[sector-1].phiMaxVec[region] );
      //
      if(m_addSectors && isNSW) {
        extendSeed(detectorHoughTransforms, truthHits, foundTruthHits, road, houghDataPerSectorVec->vec[sectorN-1] );
        associatePhiMaxima( road, houghDataPerSectorVec->vec[sectorN-1].phiMaxVec[region] );
        extendSeed(detectorHoughTransforms, truthHits, foundTruthHits, road, houghDataPerSectorVec->vec[sectorP-1] );
        associatePhiMaxima( road, houghDataPerSectorVec->vec[sectorP-1].phiMaxVec[region] );
      }

      if( road.neighbouringRegion != MuonStationIndex::DetectorRegionUnknown ) {
        associatePhiMaxima( road, houghDataPerSectorVec->vec[sector-1].phiMaxVec[road.neighbouringRegion] );
      }
      // if close to a sector boundary, try adding maxima in that sector as well
      if( road.neighbouringSector != -1 ) {
        ATH_MSG_DEBUG("  Adding neighbouring sector " << road.neighbouringSector );
        extendSeed(detectorHoughTransforms, truthHits, foundTruthHits, road, houghDataPerSectorVec->vec[road.neighbouringSector-1] );
        associatePhiMaxima( road, houghDataPerSectorVec->vec[road.neighbouringSector-1].phiMaxVec[region] );
      }

      // finally deal with the case that we have both neighbouring region and sector
      if( road.neighbouringRegion != MuonStationIndex::DetectorRegionUnknown && road.neighbouringSector != -1) {
        associatePhiMaxima( road, houghDataPerSectorVec->vec[ road.neighbouringSector-1].phiMaxVec[road.neighbouringRegion] );
      }

      // merge phi maxima
      mergePhiMaxima( road );

      // add maxima to seed exclusion list
      associatedMaxima.insert(road.maxima.begin(),road.maxima.end());
      
      ATH_MSG_DEBUG(" New road " << road.maxima.size() );
      for( auto max : road.maxima ){
	MuonStationIndex::ChIndex             chIndex = max->hough->m_descriptor.chIndex;
	MuonStationIndex::LayerIndex          layer   = Muon::MuonStationIndex::toLayerIndex(chIndex);
	MuonStationIndex::DetectorRegionIndex region  = max->hough->m_descriptor.region;
	ATH_MSG_DEBUG(" Sector " << max->hough->m_descriptor.sector << " " << Muon::MuonStationIndex::regionName(region) 
		      << " " << Muon::MuonStationIndex::layerName(layer)
		      << " maximum " << max->max << " position " << max->pos << " angle " << max->theta << " ptr " << max);
      }
      
      bool insert = true;
      for( auto& oldRoad : roads ){
        std::vector< const MuonHough::MuonLayerHough::Maximum* > intersection;
        std::set_intersection(oldRoad.maximumSet.begin(),oldRoad.maximumSet.end(),
                              road.maximumSet.begin(),road.maximumSet.end(),
                              std::back_inserter(intersection));
        unsigned int intersectionSize = intersection.size();
        unsigned int oldRoadSize      = oldRoad.maximumSet.size();
        unsigned int roadSize         = road.maximumSet.size();
        ATH_MSG_VERBOSE(" Overlap check " << intersectionSize << " old " << oldRoadSize << " new " << roadSize << " old ptr " << oldRoad.seed );
        if( intersectionSize == 0 ) continue; 
        if( intersectionSize == roadSize ) {
          insert = false; // discard
          break;
        }else if( intersectionSize == oldRoadSize ){
          oldRoad = road; // replace 
          insert = false;
          break; 
        }
      }
      
      // add road to list
      if( insert ) roads.push_back(road);
    }
    
    
  }


  void MuonLayerHoughTool::mergePhiMaxima( MuonLayerHoughTool::Road& road ) const {
    // input -> list of phiMaxima on road
    // returns some mergedPhiMaxima -> is this "summed" over layers?
    
    auto maximaSortingLambda = [road]( const MuonHough::MuonPhiLayerHough::Maximum* m1, const MuonHough::MuonPhiLayerHough::Maximum* m2 ) { 
                                   if( m1->max == m2->max ){
                                     if (m1->sector == m2->sector){   // prefer the same sector as the seed sector
                                       if( m1->hits.size() == m2->hits.size() ) {
                                         if( m1->pos == m2->pos ) {
                                           if( std::abs(m1->binposmax - m1->binposmin) == std::abs(m2->binposmax - m2->binposmin) ) {
                                             return (m1->binposmin) < (m2->binposmin);
                                           }
                                           return std::abs(m1->binposmax - m1->binposmin) < std::abs(m2->binposmax - m2->binposmin);
                                         }
                                         return m1->pos < m2->pos;
                                       }
                                       return m1->hits.size() < m2->hits.size();  // least hits -> most collimated maximum
                                     }
                                     return m1->sector == road.seed->hough->m_descriptor.sector;
                                   }
                                   return m1->max > m2->max; 
                                 };
        
    std::stable_sort(road.phiMaxima.begin(),road.phiMaxima.end(), maximaSortingLambda);
    
    ATH_MSG_VERBOSE("Merging phi maxima " << road.phiMaxima.size() );
    std::set<MuonHough::MuonPhiLayerHough::Maximum*> associatedPhiMaxima;
    for( auto pit = road.phiMaxima.begin();pit!=road.phiMaxima.end();++pit ){  // loop over phi maxima
      if( associatedPhiMaxima.count(*pit) ) continue; //check if maximum is already in associatedPhiMaxima
      associatedPhiMaxima.insert(*pit);
      MuonHough::MuonPhiLayerHough::Maximum phiMaximum = **pit;
      ATH_MSG_VERBOSE("  phi maxima " << phiMaximum.pos << " val " << phiMaximum.max );

      bool wasExtended = false;
      for( auto pit1 = pit+1;pit1!=road.phiMaxima.end();++pit1 ){
        if( (*pit1)->binposmax >= phiMaximum.binposmin && (*pit1)->binposmin <= phiMaximum.binposmax ){
          ATH_MSG_VERBOSE("    merging maxima " << phiMaximum.pos << " val " << phiMaximum.max << " " << (*pit1)->pos << " val " << (*pit1)->max );
          phiMaximum.hits.insert(phiMaximum.hits.end(),(*pit1)->hits.begin(),(*pit1)->hits.end());
          associatedPhiMaxima.insert(*pit1);
          wasExtended = true;
        }
      }

      if( wasExtended ) {
        // refind maximum
        MuonHough::MuonPhiLayerHough localHough(60, -M_PI, M_PI, ( (*pit)->hough ? (*pit)->hough->m_region : MuonStationIndex::DetectorRegionUnknown ) );
        std::vector<MuonHough::PhiHit*> hits = phiMaximum.hits;
        /* too ambiguous producing irreproducibilities because of sorting by pointer value
        std::stable_sort(hits.begin(),hits.end(),[]( const MuonHough::PhiHit* h1,
                                                     const MuonHough::PhiHit* h2 ){ return h1->layer < h2->layer; } );
        */

        std::stable_sort(hits.begin(),hits.end(),[]( const MuonHough::PhiHit* h1,
                                                     const MuonHough::PhiHit* h2 ){
                                                      if( h1->layer == h2->layer ) {
                                                        if( h1->w == h2->w ) {
                                                          if( h1->r == h2->r ) {
                                                            if( std::abs(h1->phimax - h1->phimin) == std::abs(h2->phimax - h2->phimin) ){
                                                              if( h1->phimin == h2->phimin ) return h1->phimax < h2->phimax;
                                                              return h1->phimin < h2->phimin;
                                                            }
                                                            return std::abs(h1->phimax - h1->phimin) < std::abs(h2->phimax - h2->phimin);
                                                          }
                                                          return h1->r < h2->r;
                                                        }
                                                        return h1->w > h2->w;
                                                      }
                                                      return h1->layer < h2->layer;
                                                     } );

        ATH_MSG_VERBOSE("  updating phi maximum " << phiMaximum.pos  << " bin " << phiMaximum.binpos 
                        << " val " << phiMaximum.max << " number of hits " << hits.size() );
        if( msgLvl(MSG::VERBOSE) ) localHough.setDebug(true);
        localHough.fillLayer2(hits);
        localHough.findMaximum(phiMaximum,0.9);
        localHough.associateHitsToMaximum(phiMaximum,hits);
        ATH_MSG_VERBOSE("  updated phi maxima " << phiMaximum.pos << " bin " << phiMaximum.binpos 
                        << " val " << phiMaximum.max << " number of hits " << phiMaximum.hits.size() );
        phiMaximum.hough = (*pit)->hough; // set back pointer to transform
      }
      road.mergedPhiMaxima.push_back(phiMaximum);
    }
  }

  // maximum in middle layer
  // says look in other layers
  // if yes, combine them
  // gets on road
  // roads are combinations of maxima
  
  
  void MuonLayerHoughTool::extendSeed(MuonHough::MuonDetectorHough& detectorHoughTransforms, 
      std::set<Identifier>& truthHits, std::set<Identifier>& foundTruthHits, 
      MuonLayerHoughTool::Road& road, MuonLayerHoughTool::HoughDataPerSector& sectorData ) const { //const {
    if( !road.seed ) return;
    
    RegionMaximumVec& maxVec = sectorData.maxVec;

    // gather locality information on seed
    MuonHough::MuonLayerHough::Maximum&   seed = *road.seed;
    MuonStationIndex::LayerIndex          seedLayer   = Muon::MuonStationIndex::toLayerIndex(seed.hough->m_descriptor.chIndex);
    MuonStationIndex::DetectorRegionIndex region  = seed.hough->m_descriptor.region;

    // loop over layers in the same region as the seed ( inner, middle, outer)
    for( int lay=0;lay<Muon::MuonStationIndex::LayerIndexMax;++lay ){
      MuonStationIndex::LayerIndex layer   = static_cast<MuonStationIndex::LayerIndex>(lay);
      if( layer == seedLayer && seed.hough->m_descriptor.sector == sectorData.sector ) continue; 
      
      // untrue -> look in neighboring layer
      // true -> look only in this layer
      double distanceCut = layer == seedLayer ? 500. : (double)m_extrapolationDistance;

      unsigned int layerHash = MuonStationIndex::sectorLayerHash(region,layer);
      
      // fetching vector of maxima for given region and layer
      const MaximumVec& maxima = maxVec[layerHash];
      if( maxima.empty() ) continue;

      ATH_MSG_DEBUG("Associating maxima in " << MuonStationIndex::regionName(region) << " " << MuonStationIndex::layerName(layer)
                    << " size " << maxima.size() ); 
      // loop over maxima in layer
      for( auto mit = maxima.begin();mit!=maxima.end();++mit ){
        MuonHough::MuonLayerHough::Maximum* candMaximum = *mit;
        // extrapolate seed to layer assuming a pointing straight line or parabolic
        // add maximum to road if close enough
        float yloc_diff =  MuonHough::extrapolate(seed, *candMaximum, m_doParabolicExtrapolation);
        if( fabs( MuonHough::extrapolate(seed, *candMaximum, m_doParabolicExtrapolation) ) < distanceCut ) {
          ATH_MSG_VERBOSE(" Adding maximum position " << candMaximum->pos << " intersect diff" << yloc_diff );
          road.add(candMaximum);
        }else{
          ATH_MSG_VERBOSE(" Maximum position: y " << candMaximum->pos 
                                        <<  " x " << candMaximum->hough->m_descriptor.referencePosition 
                                    << " seed y " << seed.hough->m_descriptor.referencePosition 
                                         << " x " << seed.pos << " intersect diff " << yloc_diff );
        }
      }
    }

    // check if the maximum is close to the detector boundary, if yes look for maxima in the neighbouring region, skip BarrelExtended
    if( seedLayer == MuonStationIndex::BarrelExtended ) return;

    ATH_MSG_DEBUG("Checking Barrel/Endcap overlaps: min dist edge " << seed.pos-seed.hough->m_descriptor.yMinRange 
                                               << " max dist edge " << seed.pos-seed.hough->m_descriptor.yMaxRange
                                                         << " pos " << seed.pos 
                                                       << " range " << seed.hough->m_descriptor.yMinRange << " " << seed.hough->m_descriptor.yMaxRange );

    if( fabs(seed.pos-seed.hough->m_descriptor.yMinRange) < 4000. || fabs(seed.pos-seed.hough->m_descriptor.yMaxRange) < 4000. ){

      // asumes region is barrel and looks in adjacent regions (clever logic TM here)
      MuonStationIndex::DetectorRegionIndex neighbourRegion = MuonStationIndex::Barrel;
      if( region == MuonStationIndex::Barrel ){
        if( seed.pos < 0 ) neighbourRegion = MuonStationIndex::EndcapC;
        else               neighbourRegion = MuonStationIndex::EndcapA;
      } // in all other cases the neigbourRegion is definitely barrel
      
      // looping over all layers in neigbouring region
      for( int lay=0;lay<Muon::MuonStationIndex::LayerIndexMax;++lay ){
        MuonStationIndex::LayerIndex layer   = static_cast<MuonStationIndex::LayerIndex>(lay);
        
        // skip barrel combinations with BEE
        if( region == MuonStationIndex::Barrel && layer == MuonStationIndex::BarrelExtended ) continue;

        double distanceCut = 1000.;

        // get maxima from neigboring region
        unsigned int layerHash = MuonStationIndex::sectorLayerHash(neighbourRegion,layer);
        const MaximumVec& maxima = maxVec[layerHash];
        if( maxima.empty() ) continue;
        ATH_MSG_DEBUG("Associating maxima in neighbouring region " << MuonStationIndex::regionName(neighbourRegion) 
                                                            << " " << MuonStationIndex::layerName(layer) 
                                                       << " hash " << layerHash 
                                                       << " size " << maxima.size() ); 
  
        // loop over maxima per layer
        for( auto mit = maxima.begin();mit!=maxima.end();++mit ){
          MuonHough::MuonLayerHough::Maximum* candMaximum = *mit;
          // extrapolate seed to layer assuming a pointing straight line, swap coordinates
          float yloc_diff =  MuonHough::extrapolate(seed, *candMaximum, m_doParabolicExtrapolation);
          ATH_MSG_VERBOSE(" Maximum position: y " << candMaximum->pos 
                                         << " x " << candMaximum->hough->m_descriptor.referencePosition 
                                    << " seed y " << seed.hough->m_descriptor.referencePosition 
                                         << " x " << seed.pos 
                                 << " intersect diff " << yloc_diff );

          if( fabs(yloc_diff) < distanceCut ) {
            road.add(candMaximum);
            road.neighbouringRegion = neighbourRegion;
          }
        }      
      }
    }    

    // search for phiMaxima using the etaMaximum of the road in the current sector
    std::set<const TgcClusterObj3D*> tgcClusters;
    std::set<Identifier> triggerLayers;
    const MaximumVec& maxima = road.maxima;
    for( auto mit = maxima.begin();mit!=maxima.end();++mit ){
      MuonHough::MuonLayerHough::Maximum* maximum = *mit;
      if (maximum->hough->m_descriptor.sector != sectorData.sector) continue; // skip cases where a maximum on the road does not belong to the currently examined sector

      // gather tgcClusters associated to the hits of the maxima
      for( auto ehit = maximum->hits.begin();ehit!=maximum->hits.end();++ehit ){
        const MuonHough::Hit& etaHit = **ehit;
        if (etaHit.tgc){
          if (!etaHit.tgc->phiCluster.hitList.empty()) tgcClusters.insert(etaHit.tgc);
        }
        else if( etaHit.prd ){
          triggerLayers.insert(m_idHelperSvc->gasGapId(etaHit.prd->identify()));
        }
      }
    }
      
    MuonHough::MuonPhiLayerHough& phiHough = detectorHoughTransforms.phiHough( region ); // get phi transform in the same region as the seed
            
    // gather phiHits in sector that match the etaHits of the maximum
    PhiHitVec phiHitsInMaximum;
    PhiHitVec& phiHits = sectorData.phiHitVec[region];      
    for( PhiHitVec::iterator phit = phiHits.begin(); phit != phiHits.end();++phit ){
      MuonHough::PhiHit *phiHit = *phit;
      if (phiHit->tgc){
        if (tgcClusters.find((*phit)->tgc) != tgcClusters.end()) phiHitsInMaximum.push_back(phiHit);
      }
      else if (phiHit->prd){
        if (triggerLayers.find(m_idHelperSvc->gasGapId(phiHit->prd->identify())) != triggerLayers.end()) phiHitsInMaximum.push_back(phiHit);
      }
    }
    
    // fill phi hits 
    ATH_MSG_DEBUG("extendSeed: Filling s"  << sectorData.sector
                               <<  " "  << MuonStationIndex::regionName(region)
                               << " phiHitsInMaxima " << phiHitsInMaximum.size()
                               << " phi hits:  "      << phiHits.size() );

    if( !findMaxima(truthHits, foundTruthHits, phiHough,phiHitsInMaximum,
        sectorData.phiMaxVec[region], sectorData.sector) || sectorData.phiMaxVec[region].empty() ) {
      ATH_MSG_DEBUG("extendSeed: No phi maxima found in  s" << sectorData.sector << " " << MuonStationIndex::regionName(region) );
      return;
    }

    ++sectorData.nphilayersWithMaxima[region];
    sectorData.nphimaxHitsInRegion[region] += sectorData.phiMaxVec[region].front()->max;

    ATH_MSG_DEBUG("extendSeed: Sector phiMaxima Summary:  s" << sectorData.sector <<  " " << MuonStationIndex::regionName(region) <<   " " << sectorData.nphilayersWithMaxima[region] 
                                                             <<  " -> nPhiMaxima: " << sectorData.phiMaxVec[region].size()
                                                             <<  " max sum: "    << sectorData.nphimaxHitsInRegion[region] );


  }
  
  // phi hits are not separated into inner middle outer
  // maxima found in road
  void MuonLayerHoughTool::associatePhiMaxima( MuonLayerHoughTool::Road& road, MuonLayerHoughTool::PhiMaximumVec& phiMaxima ) const {
    ATH_MSG_DEBUG("associateMaximaToPhiMaxima: phi maxima " << phiMaxima.size() );
    if( !road.seed ) return;

    // loop over phi maxima
    for( auto pit = phiMaxima.begin();pit!=phiMaxima.end();++pit ){
      
      // reference to phi maximum
      MuonHough::MuonPhiLayerHough::Maximum& pmaximum = **pit;
      
      ATH_MSG_DEBUG(" new phi maximum " << pmaximum.max << " hits " << pmaximum.hits.size() );

      // precalculate the layers + TGC clusters and add them to a set for easy access
      std::map<Identifier, std::pair<double,double> > triggerLayersPhiMinMax;
      std::map<MuonStationIndex::StIndex,std::set<const TgcClusterObj3D*> > tgcClusters;

      // loop over hits
      PhiHitVec::const_iterator phit = pmaximum.hits.begin();
      PhiHitVec::const_iterator phit_end = pmaximum.hits.end();
      for( ;phit!=phit_end;++phit ){
        const MuonHough::PhiHit& phiHit = **phit; 
        // two cases 
        // case 1: phiHit measured in TGC -> get phiHits from phiCluster
        // case 2: phiHit is prepared raw data -> use phiHit to extend the triggerLayersPhinMinMax map
        if( phiHit.tgc ){
          if(  phiHit.tgc->phiCluster.hitList.empty() ) ATH_MSG_WARNING(" TGC 3D cluster without phi hits ");
          else tgcClusters[m_idHelperSvc->stationIndex( phiHit.tgc->phiCluster.hitList.front()->identify() )].insert(phiHit.tgc);
        }
        else if( phiHit.prd ){
          Identifier gpId = m_idHelperSvc->gasGapId(phiHit.prd->identify());
          auto mit = triggerLayersPhiMinMax.find(gpId);
          if( mit == triggerLayersPhiMinMax.end() ) triggerLayersPhiMinMax[gpId] = std::make_pair(phiHit.phimin,phiHit.phimax);
          else{
            if( phiHit.phimin < mit->second.first  ) mit->second.first  = phiHit.phimin;
            if( phiHit.phimax > mit->second.second ) mit->second.second = phiHit.phimax;
          }
        }
      }
      // print out information on the triggerLayersPhiMinMax
      if( msgLvl(MSG::DEBUG) && false ){
        ATH_MSG_DEBUG("Trigger layers " << triggerLayersPhiMinMax.size() << " tgc layers " << tgcClusters.size() );
        for( auto tgcit = triggerLayersPhiMinMax.begin() ;tgcit!=triggerLayersPhiMinMax.end();++tgcit ){
          ATH_MSG_VERBOSE("  " << m_idHelperSvc->toString(tgcit->first) );
        }
        
        // loop over the stations and the contained tgcClusters found in the previous step, print out information
        std::map<MuonStationIndex::StIndex,std::set<const TgcClusterObj3D*> >::const_iterator stit = tgcClusters.begin();
        std::map<MuonStationIndex::StIndex,std::set<const TgcClusterObj3D*> >::const_iterator stit_end = tgcClusters.end();
        for( ;stit != stit_end; ++stit ){
          std::set<const TgcClusterObj3D*>::const_iterator ttit = stit->second.begin();
          std::set<const TgcClusterObj3D*>::const_iterator ttit_end = stit->second.end();
          for( ;ttit!=ttit_end;++ttit ){
            ATH_MSG_VERBOSE("  " << m_idHelperSvc->toString( (*ttit)->phiCluster.hitList.front()->identify() ) << "  nhits " <<  (*ttit)->phiCluster.hitList.size() );
          }
        }
      }

      // check if there are matching maxima in neighbouring sectors, add maximum values if confirmation is found
      // overlap counters
      int noverlaps = 0;
      int nNoOverlaps = 0;
      double phimin =  10;
      double phimax = -10;
      
      // loop over all maxima found on road
      for( auto max : road.maxima ){

        // get station information for maximum on road
        const MuonHough::MuonLayerHough::Maximum& maximum = *max;
        MuonStationIndex::StIndex stIndex = MuonStationIndex::toStationIndex(maximum.hough->m_descriptor.chIndex);

        // loop over eta hits
        for( auto ehit = maximum.hits.begin();ehit!=maximum.hits.end();++ehit ){
          const MuonHough::Hit& etaHit = **ehit; 
          if( etaHit.tgc ) {
            if( etaHit.tgc->etaCluster.hitList.empty() ) ATH_MSG_WARNING(" TGC 3D cluster without eta hits " );
            else{
              if( tgcClusters[stIndex].count(etaHit.tgc) ){
                // now loop over phi maximum and find phi hit
                for( auto phit = pmaximum.hits.begin(); phit!=pmaximum.hits.end();++phit ){
                  const MuonHough::PhiHit& phiHit = **phit; 
                  if( phiHit.tgc == etaHit.tgc ){
                    if( phiHit.phimin < phimin ) phimin = phiHit.phimin;
                    if( phiHit.phimin > phimax ) phimax = phiHit.phimax;
                    break;
                  }
                }
                ++noverlaps;
              }
              else{
                ++nNoOverlaps;
              }
            }
          }else if( etaHit.prd ){
            if( !m_idHelperSvc->isRpc(etaHit.prd->identify()) && !m_idHelperSvc->issTgc(etaHit.prd->identify())) continue;
            Identifier gpId = m_idHelperSvc->gasGapId( etaHit.prd->identify() );
            auto mit = triggerLayersPhiMinMax.find(gpId);
            if( mit == triggerLayersPhiMinMax.end() )  ++nNoOverlaps;
            else{
              if( mit->second.first  < phimin ) phimin = mit->second.first;
              if( mit->second.second > phimax ) phimax = mit->second.second;
              ++noverlaps;
            }
          }
        } // loop over hits in maximum
      } // loop over maxima in road
      
      // if overlaps are found, add the phi maximum in question to the road
      if( noverlaps > 0 ) {
        road.add(&pmaximum);
        // check if we are close to a sector boundary
        std::vector<int> sectors;
        m_sectorMapping.getSectors(phimin,sectors);
        if( sectors.size() > 1 ){
          for( auto sec : sectors ){
            if( sec != road.seed->hough->m_descriptor.sector ) road.neighbouringSector = sec;
          }
        }else{
          std::vector<int> sectors;
          m_sectorMapping.getSectors(phimax,sectors);
          if( sectors.size() > 1 ){
            for( auto sec : sectors ){
              if( sec != road.seed->hough->m_descriptor.sector ) road.neighbouringSector = sec;
            }
          }
        }
      }
      ATH_MSG_DEBUG(" Overlap with Phi maximum: overlap " << noverlaps << " no overlap " << nNoOverlaps 
                    << " phimin " << phimin << " phimax " << phimax << " neighbouring sector " << road.neighbouringSector );  

    }
  }

  // takes the maxima from a given sector and tries to associate it with the maxima of the adjacent sectors
  void MuonLayerHoughTool::associateMaximaInNeighbouringSectors( MuonLayerHoughTool::HoughDataPerSector& houghData, std::vector<MuonLayerHoughTool::HoughDataPerSector>& houghDataPerSectorVec ) const {
    
    ATH_MSG_DEBUG(" looping over eta maxima");

    // now loop over eta maxima per layer
    for( unsigned int regLay=0;regLay<houghData.maxVec.size();++regLay ){
  
      MaximumVec& maxima = houghData.maxVec[regLay];
      int sector = houghData.sector;

      // loop over two neighbouring sectors
      for( int i=0;i<2;++i ){
  
        // calculate neighbouring sector index
        int sectorN = (i == 0) ? sector - 1 : sector + 1;
        if( i == 0 && sector == 1 ) sectorN = 16;
        if( i == 1 && sector == 16 ) sectorN = 1;
        
        MuonLayerHoughTool::HoughDataPerSector& houghDataN = houghDataPerSectorVec[sectorN-1];
        
        MaximumVec& maximaN = houghDataN.maxVec[regLay];
        
        // loop over maxima in layer
        MaximumVec::iterator mit     = maxima.begin();
        MaximumVec::iterator mit_end = maxima.end();
        for( ;mit!=mit_end;++mit ){
        
          // reference to maximum
          MuonHough::MuonLayerHough::Maximum& maximum = **mit;
          
          if( !maximum.hough ){
            ATH_MSG_WARNING("Maximum without associated hough transform! " );
            continue;
          }
        
          // loop over maxima per layer in neighbouring sector
          MaximumVec::iterator nmit = maximaN.begin();
          MaximumVec::iterator nmit_end = maximaN.end();
          for( ;nmit!=nmit_end;++nmit ){
        
            // reference to maximum
            MuonHough::MuonLayerHough::Maximum& maximumN = **nmit;
            if( !maximumN.hough ){
              ATH_MSG_WARNING("Maximum without associated hough transform! " );
              continue;
            }
            
            // calculate the position of the first maximum in the reference frame of the other sector
            double rcor = maximumN.hough->m_descriptor.referencePosition*m_sectorMapping.transformRToNeighboringSector( maximum.pos,sector,sectorN)/maximum.hough->m_descriptor.referencePosition;
            double dist = rcor - maximumN.pos;
            ATH_MSG_DEBUG(" maximumN.hough " << maximumN.hough->m_descriptor.referencePosition << " maximum.hough " << maximum.hough->m_descriptor.referencePosition << " maximumN.pos " << maximumN.pos << " maximum.pos " << maximum.pos << rcor << " distance " << dist );
            if( fabs(dist) > 100 ) continue;
            houghData.maxAssociationMap[&maximum].push_back(&maximumN);
            houghDataN.associatedToOtherSector.insert(&maximumN);
        
            ATH_MSG_DEBUG(" Found maximum in neighbouring sector: max " << maximum.max
              << " pos " << rcor << " maxN " <<  maximumN.max << " pos " << maximumN.pos 
              << " distance " << dist   );
            
            // loop over first and second maximum
            for( int nn = 0; nn < 2; ++nn ){
        
              // setting info for the debug-info objects of the hits
              MuonHough::MuonLayerHough::Maximum& maxi = nn == 0 ? maximum : maximumN;
              MuonHough::MuonLayerHough::Maximum& maxi2 = nn == 0 ? maximumN : maximum;
              ATH_MSG_VERBOSE(" Maximum " << nn << " hits " << maxi.hits.size() );
              HitVec::const_iterator ehit = maxi.hits.begin();
              HitVec::const_iterator ehit_end = maxi.hits.end();
              for( ;ehit!=ehit_end;++ehit ){
                MuonHough::Hit& hit = **ehit; 
                if( hit.debugInfo() ) {
                  hit.debugInfo()->phn = maxi2.max;
                  Identifier id = hit.tgc ? hit.tgc->etaCluster.hitList.front()->identify() : hit.prd->identify();
                  ATH_MSG_VERBOSE(" " << m_idHelperSvc->toString(id) << " setphn " << hit.debugInfo()->phn);
                }
              }
            }
          }
        }
      } 
    }
  }


  void MuonLayerHoughTool::associateMaximaToPhiMaxima( MuonStationIndex::DetectorRegionIndex region, MuonLayerHoughTool::HoughDataPerSector& houghData,
         std::map< MuonHough::MuonPhiLayerHough::Maximum*, MuonLayerHoughTool::MaximumVec >& phiEtaAssociations,
         MuonLayerHoughTool::RegionMaximumVec& unassEtaMaxima) const
  {
         
    std::set<MuonHough::MuonLayerHough::Maximum*> associatedMaxima;

    PhiMaximumVec& phiMaxima = houghData.phiMaxVec[region];

    ATH_MSG_DEBUG("associateMaximaToPhiMaxima in sector " << houghData.sector << ": phi maxima " << phiMaxima.size() ); // !!!!
    // loop over phi maxima
    PhiMaximumVec::iterator pit = phiMaxima.begin(); 
    PhiMaximumVec::iterator pit_end = phiMaxima.end(); 
    for( ;pit!=pit_end;++pit ){
      
      // reference to phi maximum
      MuonHough::MuonPhiLayerHough::Maximum& phiMaximum = **pit;
      
      ATH_MSG_DEBUG(" Considering phi maximum " << phiMaximum.max << " hits " << phiMaximum.hits.size() ); // !!!!

      // store associated maxima
      MaximumVec associatedMaximaVec; // !!!!

      // precalculate the layers + TGC clusters and add them to a set for easy access
      // std::map< Identifier,std::pair<double,double> > triggerLayersP;
      std::set<Identifier> triggerLayers;
      std::map<MuonStationIndex::StIndex,std::set<const TgcClusterObj3D*> > tgcClusters;
      
      // loop over hits
      PhiHitVec::const_iterator phit = phiMaximum.hits.begin();
      PhiHitVec::const_iterator phit_end = phiMaximum.hits.end();
      for( ;phit!=phit_end;++phit ){
        const MuonHough::PhiHit& phiHit = **phit;
              // double rmin = 0.;
              // double rmax = 0.;
        
        if( phiHit.tgc ){
          if(  phiHit.tgc->phiCluster.hitList.empty() ) ATH_MSG_WARNING(" TGC 3D cluster without phi hits ");
          else tgcClusters[m_idHelperSvc->stationIndex( phiHit.tgc->phiCluster.hitList.front()->identify() )].insert(phiHit.tgc);
        }
        else if( phiHit.prd ){
          Identifier colId = phiHit.prd->identify();
          Identifier layId = m_idHelperSvc->gasGapId( colId );
          triggerLayers.insert(layId);
        }
            }
            if( msgLvl(MSG::DEBUG) ){
              ATH_MSG_DEBUG("Trigger layers " << triggerLayers.size() << " tgc layers " << tgcClusters.size() );
              auto tgcit = triggerLayers.begin();
              auto tgcit_end = triggerLayers.end();
              for( ;tgcit!=tgcit_end;++tgcit ){
                ATH_MSG_VERBOSE("  " << m_idHelperSvc->toString(*tgcit) );
              }
              
              std::map<MuonStationIndex::StIndex,std::set<const TgcClusterObj3D*> >::const_iterator stit = tgcClusters.begin();
              std::map<MuonStationIndex::StIndex,std::set<const TgcClusterObj3D*> >::const_iterator stit_end = tgcClusters.end();
              for( ;stit != stit_end; ++stit ){
                std::set<const TgcClusterObj3D*>::const_iterator ttit = stit->second.begin();
                std::set<const TgcClusterObj3D*>::const_iterator ttit_end = stit->second.end();
                for( ;ttit!=ttit_end;++ttit ){
                  ATH_MSG_VERBOSE("  " << m_idHelperSvc->toString( (*ttit)->phiCluster.hitList.front()->identify() ) << "  nhits " <<  (*ttit)->phiCluster.hitList.size() );
                }
              }
            }
        
            ATH_MSG_DEBUG(" looping over eta maxima");
        
            // now loop over eta maxima per layer
        for( unsigned int lay=0;lay<MuonStationIndex::LayerIndexMax;++lay ){
          MuonStationIndex::LayerIndex layer = static_cast<MuonStationIndex::LayerIndex>(lay);
          unsigned int layerHash = MuonStationIndex::sectorLayerHash(region,layer);
          MaximumVec& maxima = houghData.maxVec[layerHash];
          if( maxima.empty() ) continue;
          MuonStationIndex::StIndex stIndex = MuonStationIndex::toStationIndex(region,layer);
        
        // loop over maxima per layer
        MaximumVec::iterator mit = maxima.begin();
        MaximumVec::iterator mit_end = maxima.end();
        for( ;mit!=mit_end;++mit ){
        
          // skip maxima that were already associated to a neighbouring sector
          if( houghData.associatedToOtherSector.count(*mit) ) continue;
        
          // reference to maximum
          MuonHough::MuonLayerHough::Maximum& maximum = **mit;
          
          // check if there are matching maxima in neighbouring sectors, add maximum values if confirmation is found
          int totmax = 0;
          int ntrigconfirm = 0;
          MaximumAssociationMap::iterator pos = houghData.maxAssociationMap.find(&maximum);
          if( pos != houghData.maxAssociationMap.end() ){
            for( MaximumVec::iterator amit = pos->second.begin();amit !=pos->second.end();++amit ){
              if( (*amit)->max > totmax ) totmax = (*amit)->max;
              ntrigconfirm += (*amit)->triggerConfirmed;
            }
          }
          totmax += maximum.max;
          ntrigconfirm += maximum.triggerConfirmed;
        
          ATH_MSG_DEBUG("   new eta maximum " << MuonStationIndex::stName(stIndex) << " val " << maximum.max 
            << " neightbour confirmed value " << totmax << " trigger confirmations " << ntrigconfirm );
          
          
          // overlap counters
          int nmmHits = 0;
          int ntgcOverlaps = 0;
          int nrpcOverlaps = 0;
          int nstgcOverlaps = 0;
          int ntgcNoOverlaps = 0;
          int nrpcNoOverlaps = 0;
          int nstgcNoOverlaps = 0;
        
          // loop over hits
          HitVec::const_iterator ehit = maximum.hits.begin();
          HitVec::const_iterator ehit_end = maximum.hits.end();
          for( ;ehit!=ehit_end;++ehit ){
            const MuonHough::Hit& etaHit = **ehit; 
            if( etaHit.tgc ) {
              if( etaHit.tgc->etaCluster.hitList.empty() ) ATH_MSG_WARNING(" TGC 3D cluster without eta hits " );
              else{
                if( tgcClusters[stIndex].count(etaHit.tgc) ) ++ntgcOverlaps;
                else                                         ++ntgcNoOverlaps;
              }
            }
            else if( etaHit.prd ){
              Identifier layId = m_idHelperSvc->gasGapId( etaHit.prd->identify() );
              ATH_MSG_VERBOSE(" eta layer hit " << m_idHelperSvc->toString(layId) );
              if( m_idHelperSvc->isMM(layId) ) ++nmmHits;
              if( triggerLayers.count(layId) ){
                if( m_idHelperSvc->isRpc(layId) )       ++nrpcOverlaps;
                else if( m_idHelperSvc->issTgc(layId) ) ++nstgcOverlaps;
            }
            else{
              if( m_idHelperSvc->isRpc(layId) )       ++nrpcNoOverlaps;
              else if( m_idHelperSvc->issTgc(layId) ) ++nstgcNoOverlaps;
            }
          }
        }
    
      // cuts on NSW endcap only for now
      if( nmmHits + nstgcNoOverlaps + nstgcOverlaps > 0 ) {
      // select 
      if( maximum.pos < 1200. ){
        if( totmax < 8 )  {
          ATH_MSG_DEBUG("  maximum failed cut " << totmax  << " cut 8, position " << maximum.pos );
          continue;
        }
      }
      else if( maximum.pos > 4300. ){
        if( totmax < 8 )  {
          ATH_MSG_DEBUG("  maximum failed cut " << totmax  << " cut 8, position " << maximum.pos );
          continue;
        }
      }
      else{
        if( totmax < 12 ) {
          ATH_MSG_DEBUG("  maximum failed cut " << totmax  << " cut 12, position " << maximum.pos );
          continue;
        }
      }
    }
     
    if( m_ntuple ) {
      m_ntuple->fill(nstgcOverlaps,nstgcNoOverlaps);
    }

    
    ATH_MSG_DEBUG(" Overlap with Phi maximum: tgc " << ntgcOverlaps << " stgc " << nstgcOverlaps << " rpc " << nrpcOverlaps
      << " nphiTgc " << tgcClusters[stIndex].size() << " trigLay " << triggerLayers.size() );  
    if( stIndex == MuonStationIndex::EM && !tgcClusters[stIndex].empty() && ntgcOverlaps == 0 ) {
      ATH_MSG_VERBOSE(" No association in StationLayer " << MuonStationIndex::stName(stIndex) << " tgcs overlaps " << ntgcOverlaps 
          << " on phi maximum " << tgcClusters[stIndex].size() );  
      continue;
    }
    if( stIndex == MuonStationIndex::EI && !tgcClusters[stIndex].empty() && ntgcOverlaps == 0 ) {
      ATH_MSG_VERBOSE(" No association in StationLayer " << MuonStationIndex::stName(stIndex) << " tgcs overlaps " << ntgcOverlaps 
          << " on phi maximum " << tgcClusters[stIndex].size() );  
      continue;
    }
    if( stIndex == MuonStationIndex::EI && nstgcOverlaps == 0 && nstgcNoOverlaps != 0 ) {
      ATH_MSG_VERBOSE(" No association in StationLayer " << MuonStationIndex::stName(stIndex) << " stgcs without overlaps " << nstgcNoOverlaps );  
      continue;
    }
    // require STGC confirmation 
    if( m_requireTriggerConfirmationNSW && nmmHits > 0 && ntrigconfirm == 0 ) continue;
    
    associatedMaxima.insert(&maximum);
    associatedMaximaVec.push_back(&maximum);

    // check if there are matching maxima in neighbouring sectors
    if( pos != houghData.maxAssociationMap.end() ){
      associatedMaxima.insert(pos->second.begin(),pos->second.end());
      associatedMaximaVec.insert(associatedMaximaVec.end(),pos->second.begin(),pos->second.end());
    }
  }
      }
      
      if( associatedMaximaVec.empty() ) continue;
      ATH_MSG_DEBUG(" processed phi maximum, associated eta maxima " << associatedMaximaVec.size() ); 
      phiEtaAssociations[*pit] = associatedMaximaVec;
      
    }
    
    // finally idenitify all unassociated maxima and add them to the unassociated maxima list
    // now loop over eta maxima per layer
    for( unsigned int lay=0;lay<MuonStationIndex::LayerIndexMax;++lay ){
      MuonStationIndex::LayerIndex layer = static_cast<MuonStationIndex::LayerIndex>(lay);
      unsigned int layerHash = MuonStationIndex::sectorLayerHash(region,layer);

      if( layer >= (int)unassEtaMaxima.size() ){
        ATH_MSG_WARNING(" size of unassEtaMaxima too small for region " << unassEtaMaxima.size() << " region " << MuonStationIndex::regionName(region) );
        break;
      }
      MaximumVec& maxima = houghData.maxVec[layerHash];
 
      // loop over maxima per layer
      MaximumVec::iterator mit = maxima.begin();
      MaximumVec::iterator mit_end = maxima.end();
      for( ;mit!=mit_end;++mit ){
        if( associatedMaxima.count(*mit) ) continue;
        unassEtaMaxima[layer].push_back(*mit);
        ATH_MSG_DEBUG(" unassociated maximum in layer " << MuonStationIndex::layerName(layer) << " max-val " << (*mit)->max );
      }
    }
  }



  void MuonLayerHoughTool::createPatternCombinations( MuonLayerHoughTool::RegionMaximumVec& maxima,
                  MuonPatternCombinationCollection& patternCombis )  const {
    
    ATH_MSG_DEBUG("Creating pattern combinations for eta patterns ");

    std::vector<MuonPatternChamberIntersect> chamberData;

    // bool isEndcap = maxima.size() == 5;

    // loop over layers
    RegionMaximumVec::iterator lit = maxima.begin();
    RegionMaximumVec::iterator lit_end = maxima.end();
    for( ;lit!=lit_end;++lit ){

      // create vector for prds per chamber
      std::map< Identifier, std::set< const Trk::PrepRawData* > > prdsPerChamber;

      // loop over maxima per layer
      MaximumVec::iterator mit = lit->begin();
      MaximumVec::iterator mit_end = lit->end();
      for( ;mit!=mit_end;++mit ){
  
  MuonHough::MuonLayerHough::Maximum& max = **mit;
  ATH_MSG_DEBUG("  new maximum  " << max.max << " hits " << max.hits.size() );

  // sanity check
  if( max.hits.empty() ) {
    ATH_MSG_WARNING(" Maximum without hits  ");
    continue;
  }
  ATH_MSG_DEBUG("  adding hits " << max.hits.size());

  
  // loop over hits in maximum and add them to the hit list
  HitVec::const_iterator hit = max.hits.begin();
  HitVec::const_iterator hit_end = max.hits.end();
  for( ;hit!=hit_end;++hit ) {
    Identifier chId;
    if( (*hit)->tgc ){
      chId = m_idHelperSvc->chamberId( (*hit)->tgc->etaCluster.hitList.front()->identify() );
      prdsPerChamber[chId].insert((*hit)->tgc->etaCluster.hitList.begin(),(*hit)->tgc->etaCluster.hitList.end());
    }
    else if( (*hit)->prd ){
      chId = m_idHelperSvc->chamberId( (*hit)->prd->identify() );
      prdsPerChamber[chId].insert((*hit)->prd);
    }
  }
      }

      auto sortPrdIds = [](const Trk::PrepRawData* prd1,const Trk::PrepRawData* prd2 ) { return prd1->identify() < prd2->identify(); };
      std::map< Identifier, std::set< const Trk::PrepRawData* > >::iterator chit = prdsPerChamber.begin();
      std::map< Identifier, std::set< const Trk::PrepRawData* > >::iterator chit_end = prdsPerChamber.end();
      for( ;chit!=chit_end;++chit ){
  ATH_MSG_DEBUG("Adding chamber " << m_idHelperSvc->toStringChamber(chit->first) << " hits " << chit->second.size() );
  std::vector<const Trk::PrepRawData*> prds;
  prds.insert(prds.end(),chit->second.begin(),chit->second.end());
        std::stable_sort(prds.begin(),prds.end(),sortPrdIds);
  const Trk::PrepRawData& prd = **prds.begin();
  Amg::Vector3D gpos = prd.detectorElement()->surface(prd.identify()).center();
  // create intersection and add it to combination
        ATH_MSG_DEBUG("Adding chamber with intersect phi direction " << gpos.phi() << " theta " << gpos.theta() ); 
  MuonPatternChamberIntersect intersect(gpos,gpos.unit(),prds);
  chamberData.push_back(intersect);
      }
    }
    if( chamberData.empty() ) return;

    MuonPatternCombination* combi = new MuonPatternCombination(0,chamberData);

    ATH_MSG_DEBUG(" creating new unassociated " << m_printer->print( *combi ) );
    patternCombis.push_back(combi);
  }

  void MuonLayerHoughTool::createPatternCombinations(std::set<Identifier>& truthHits, std::set<Identifier>& outputTruthHits, 
                  std::map< MuonHough::MuonPhiLayerHough::Maximum*, MuonLayerHoughTool::MaximumVec >& phiEtaAssociations,
                  MuonPatternCombinationCollection& patternCombis ) const {
    
    ATH_MSG_DEBUG("Creating pattern combinations from eta/phi combinations " << phiEtaAssociations.size() );

    // loop over the phi maxima
    std::map< MuonHough::MuonPhiLayerHough::Maximum*, MaximumVec >::const_iterator pit = phiEtaAssociations.begin();
    std::map< MuonHough::MuonPhiLayerHough::Maximum*, MaximumVec >::const_iterator pit_end = phiEtaAssociations.end();
    for( ;pit!=pit_end; ++pit ){

      if( pit->second.empty() ) continue;

      // collect phi hits per chamber
      std::map<Identifier,std::set<const Trk::PrepRawData*> > phiHitsPerChamber;
      
      // loop over hits
      PhiHitVec::const_iterator phit = pit->first->hits.begin();
      PhiHitVec::const_iterator phit_end = pit->first->hits.end();
      for( ;phit!=phit_end;++phit ){
  const MuonHough::PhiHit& hit = **phit; 
  if( hit.tgc ){
    Identifier chId = m_idHelperSvc->chamberId( hit.tgc->phiCluster.hitList.front()->identify() );
    phiHitsPerChamber[chId].insert(hit.tgc->phiCluster.hitList.begin(),hit.tgc->phiCluster.hitList.end());
  }
  else if( hit.prd ){
    Identifier chId = m_idHelperSvc->chamberId( hit.prd->identify() );
    phiHitsPerChamber[chId].insert(hit.prd);
  }
      }

      // create chamber intersections
      std::vector<MuonPatternChamberIntersect> chamberData;
      std::set<Identifier> addedPhiHits;

      // create vector for prds per chamber
      std::map< Identifier, std::set< const Trk::PrepRawData* > >      prdsPerChamber;

      // store position and direction of the first maximum in the chamber layer
      std::map< MuonStationIndex::ChIndex, std::pair<Amg::Vector3D,Amg::Vector3D> > directionsPerChamberLayer;

      // loop over eta maxima
      MaximumVec::const_iterator mit = pit->second.begin();
      MaximumVec::const_iterator mit_end = pit->second.end();
      for( ;mit!=mit_end;++mit ){
  

  const MuonHough::MuonLayerHough::Maximum& max = **mit;
  ATH_MSG_DEBUG("  new maximum  " << max.max << " hits " << max.hits.size() );

  if( !max.hough ){
    ATH_MSG_WARNING("Maximum without associated Hough Transform");
  }

  // sanity check
  if( max.hits.empty() ) {
    ATH_MSG_WARNING(" Maximum without hits  ");
    continue;
  }
  ATH_MSG_DEBUG("  adding hits " << max.hits.size());

  // loop over hits in maximum and add them to the hit list
  HitVec::const_iterator hit = max.hits.begin();
  HitVec::const_iterator hit_end = max.hits.end();
  for( ;hit!=hit_end;++hit ) {
    Identifier chId;
    if( (*hit)->tgc ){
      chId = m_idHelperSvc->chamberId( (*hit)->tgc->etaCluster.hitList.front()->identify() );
      prdsPerChamber[chId].insert((*hit)->tgc->etaCluster.hitList.begin(),(*hit)->tgc->etaCluster.hitList.end());
    }
    else if( (*hit)->prd ){
      chId = m_idHelperSvc->chamberId( (*hit)->prd->identify() );
      prdsPerChamber[chId].insert((*hit)->prd);
    }else{
      ATH_MSG_WARNING("Hit without associated PRDs");
      continue;
    }
    // the first time we have a maximun in this layer store the position and direction 
    MuonStationIndex::ChIndex chIndex = m_idHelperSvc->chamberIndex(chId);
    if( !directionsPerChamberLayer.count(chIndex) ) {
      // eta maximum has z(r) and theta parameters but these are local
      double maxpos = max.pos;
      double refPlane = 0.;
      bool isBarrel = !m_idHelperSvc->isEndcap(chId) || chIndex == MuonStationIndex::BEE;
      if( max.hough ) refPlane = max.hough->m_descriptor.referencePosition;
      else{
        if( (*hit)->tgc ) refPlane = (*hit)->tgc->p11.z();
        else{
    if( isBarrel )  refPlane = (*hit)->prd->detectorElement()->surface((*hit)->prd->identify()).center().perp();
    else            refPlane = (*hit)->prd->detectorElement()->surface((*hit)->prd->identify()).center().z();
        }              
      }
      double r = isBarrel ? refPlane : maxpos;
      double z = isBarrel ? maxpos   : refPlane;
      double theta  = max.theta;

      // go to global
      double sign = 1.;
      if(isBarrel) {
        theta += M_PI_2;
        sign = -1.;
      }

      // phi maximum has one phi from position assume global Phi definition
      double phi = pit->first->pos;//phiCor(pit->first->pos,pit->first->sector,false);
      
      CxxUtils::sincos scphi(phi);
      double sinphi = scphi.sn; 
      double cosphi = scphi.cs;

      CxxUtils::sincos sctheta(theta);
      double sintheta = sctheta.sn; 
      double costheta = sctheta.cs;

      std::pair<Amg::Vector3D,Amg::Vector3D>& posDir = directionsPerChamberLayer[chIndex];
      posDir.first  = Amg::Vector3D(r*cosphi,r*sinphi,z);
      posDir.second = Amg::Vector3D(sign*cosphi*costheta,sign*sinphi*costheta,sintheta);
      ATH_MSG_DEBUG( MuonStationIndex::chName(chIndex) 
        << " setting position: perp " << posDir.first.perp() << " z " << posDir.first.z() << " phi pos " << posDir.first.phi()
        << " direction phi  " << posDir.second.phi() << " theta pos " << posDir.first.theta()  << " direction theta " << posDir.second.theta() 
        << " ref perp " << r << " z " << z << " phi " << phi << " theta " << theta);
      if( posDir.first.dot(posDir.second) < 0. ) {
        ATH_MSG_WARNING(" direction not pointing to IP " << posDir.first.unit().dot(posDir.second));
      } 
    }

    std::map<Identifier,std::set<const Trk::PrepRawData*> >::iterator pos = phiHitsPerChamber.find(chId);
    if( pos != phiHitsPerChamber.end() ){
      std::pair<std::set<Identifier>::iterator,bool> ipos = addedPhiHits.insert(chId);
      if( ipos.second ){
        prdsPerChamber[chId].insert(pos->second.begin(),pos->second.end());
      }
    }
  }
      }

      auto sortPrdIds = [](const Trk::PrepRawData* prd1,const Trk::PrepRawData* prd2 ) { return prd1->identify() < prd2->identify(); };
      std::map< Identifier, std::set< const Trk::PrepRawData* > >::iterator chit = prdsPerChamber.begin();
      std::map< Identifier, std::set< const Trk::PrepRawData* > >::iterator chit_end = prdsPerChamber.end();
      for( ;chit!=chit_end;++chit ){
  ATH_MSG_DEBUG("Adding chamber " << m_idHelperSvc->toStringChamber(chit->first) << " hits " << chit->second.size() );
  std::vector<const Trk::PrepRawData*> prds;
  prds.insert(prds.end(),chit->second.begin(),chit->second.end());
        std::stable_sort(prds.begin(),prds.end(),sortPrdIds);
  const Trk::PrepRawData& prd = **prds.begin();

  MuonStationIndex::ChIndex chIndex = m_idHelperSvc->chamberIndex(prd.identify());
  std::map< MuonStationIndex::ChIndex, std::pair<Amg::Vector3D,Amg::Vector3D> >::const_iterator pos = directionsPerChamberLayer.find(chIndex);
  Amg::Vector3D  gpos;
  Amg::Vector3D  gdir;
  if( pos != directionsPerChamberLayer.end() ) {
    gpos = pos->second.first;
    gdir = pos->second.second;
  }else{    
    ATH_MSG_WARNING("No global position and direction found, calculating from surface");
    gpos = prd.detectorElement()->surface(prd.identify()).center();
    gdir = -1*gpos.unit();
  }

  ATH_MSG_DEBUG( "Creating intersection " <<  MuonStationIndex::chName(chIndex) 
           << " setting position: perp " << gpos.perp() << " z " << gpos.z() << " phi pos " << gpos.phi()
           << " direction phi " << gdir.phi() << " theta pos " << gpos.theta() << " theta " << gdir.theta() << " hits " << prds.size() );

  // create intersection and add it to combination
  MuonPatternChamberIntersect intersect(gpos,gdir,prds);
  chamberData.push_back(intersect);
  
  if( m_doTruth ){
    for( std::vector<const Trk::PrepRawData*>::iterator it=prds.begin();it!=prds.end();++it ){
      if( truthHits.count((*it)->identify()) ) outputTruthHits.insert((*it)->identify());
      if( !m_truthSummaryTool.empty() ) m_truthSummaryTool->add((*it)->identify(),2);
    }
  }
      }
      if( chamberData.empty() ) continue;
      if( addedPhiHits.empty() ){
  ATH_MSG_DEBUG("No phi hits selected, skipping combi ");
  continue;
      }
      MuonPatternCombination* combi = new MuonPatternCombination(0,chamberData);
      ATH_MSG_DEBUG("adding pattern combination with chambers " << chamberData.size() << " phi layers " << addedPhiHits.size()
        << std::endl << m_printer->print(*combi) );
      patternCombis.push_back(combi);
    }
  }

  bool MuonLayerHoughTool::findMaxima( std::set<Identifier>& truthHits,
                                       std::set<Identifier>& foundTruthHits,
                                       MaximumVec& seedMaxima,
                                       MuonHough::MuonLayerHough& hough,
                                       MuonLayerHoughTool::HitVec& hits, 
                                       MuonLayerHoughTool::MaximumVec& maxima ) const {
    if( hits.empty() ) return false;

    if( hough.m_descriptor.chIndex < 0 || hough.m_descriptor.chIndex >= Muon::MuonStationIndex::ChIndexMax ){
      Identifier id = hits.front()->tgc ? hits.front()->tgc->etaCluster.hitList.front()->identify() : hits.front()->prd->identify();
      ATH_MSG_WARNING("Bad ChIndex " << m_idHelperSvc->toString(id) << "  " << hough.m_descriptor.chIndex );
      return false;
    }

    // populate hough transform with hits
    std::stable_sort(hits.begin(),hits.end(),MuonHough::SortHitsPerLayer());
    if( m_debugHough ) hough.setDebug(true);
    hough.fillLayer2(hits);
      
    if( m_ntuple ) {
      updateHits(hits,hough);
    }

    
//    Muon::MuonStationIndex::StIndex stIndex = Muon::MuonStationIndex::toStationIndex(hough.m_descriptor.chIndex);
    unsigned int nmaxima = 0;
    while( nmaxima < 5 ){
      MuonHough::MuonLayerHough::Maximum maximum;
      if( hough.findMaximum( maximum, m_selectorsLoose[hough.m_descriptor.chIndex] ) ) {
        hough.associateHitsToMaximum(maximum,hits);
        ATH_MSG_DEBUG("findMaxima: Found Eta Maximum " << nmaxima 
                        <<     "  "          << maximum.max 
                        << " trigConfirmed " << maximum.triggerConfirmed 
                        << " pos "           << maximum.pos 
                        << " theta "         << maximum.theta
                        << " binPos "        << maximum.binpos 
                        << " binRange "      << maximum.binposmin << " -- " << maximum.binposmax
                        << " binTheta "      << maximum.bintheta 
                        << " nHits "         << maximum.hits.size() );

        int nmdt = 0;
        int nmm = 0;
        int ntgc = 0;
        int nstgc = 0;

        const unsigned int nHitsInMaximum = maximum.hits.size();
        for( unsigned int i = 0; i < nHitsInMaximum; ++i){
          MuonHough::Hit& hit = *(maximum.hits[i]);
          Identifier id = hit.tgc ? hit.tgc->etaCluster.hitList.front()->identify() : hit.prd->identify();
          int nhits = hit.tgc ? hit.tgc->etaCluster.hitList.size() : 1;

          if( m_idHelperSvc->isMdt(id) ) ++nmdt;
          else if( m_idHelperSvc->isTgc(id) ) ++ntgc;
          else if( m_idHelperSvc->issTgc(id) ) ++nstgc;
          else if( m_idHelperSvc->isMM(id) ) ++nmm;

          if( m_doTruth ){
            if( !m_truthSummaryTool.empty() ) m_truthSummaryTool->add(id,1);
            if( truthHits.count(id) )       foundTruthHits.insert(id);
          }

          ATH_MSG_VERBOSE("findMaxima: hit " << hit.layer << "  " << m_idHelperSvc->toString(id) << " hits " << nhits );
        }

        // only store maxima that have MDT hits        
        if( nmdt > 0 || (nmm + nstgc) > 0) {
          maxima.push_back( new MuonHough::MuonLayerHough::Maximum(maximum) );
          // add to seed list if 
          if( maximum.max > m_selectors[hough.m_descriptor.chIndex].getCutValue(maximum.pos) ) seedMaxima.push_back(maxima.back());        
          ++nmaxima;
        }
        hough.fillLayer2(maximum.hits,true);
      }
      else{
        if( nmaxima > 0 ) {
          ATH_MSG_VERBOSE("findMaxima: No more maxima found " << nmaxima );
        }
        // ?!? if nmaximo == 0 here the function should return false, I think
        break;      
      }
    }
    return true;
  }

  bool MuonLayerHoughTool::findMaxima( std::set<Identifier>& truthHits,
                                       std::set<Identifier>& foundTruthHits,
                                       MuonHough::MuonPhiLayerHough& hough,
                                       MuonLayerHoughTool::PhiHitVec& hits, 
                                       MuonLayerHoughTool::PhiMaximumVec& maxima,
                                       int sector ) const{
    if( hits.empty() ) return false;

    std::stable_sort(hits.begin(),hits.end(),MuonHough::SortHitsPerLayer());
    if( m_debugHough ) hough.setDebug(true);
    hough.fillLayer2(hits);
    
    if( m_ntuple ) {
      updateHits(hits,hough);
    }
    
    unsigned int nmaxima = 0;
    while( nmaxima < 5 ){
      MuonHough::MuonPhiLayerHough::Maximum maximum;
      if( hough.findMaximum( maximum, 1.9 ) ) {
        hough.associateHitsToMaximum(maximum,hits);

        ATH_MSG_DEBUG("findMaxima(Phi): Found Phi maximum " << nmaxima 
                                  << " height "             << maximum.max 
                                  << " pos "                << maximum.pos
                                  << " bin pos "            << maximum.binpos 
                                  << " binRange "           << maximum.binposmin << " -- " << maximum.binposmax
                                  << " nHits "              << maximum.hits.size());

        const unsigned int nHitsInMaximum = maximum.hits.size();
        for( unsigned int i = 0;i < nHitsInMaximum; ++i){
          MuonHough::PhiHit& hit = *(maximum.hits[i]);
          Identifier id = hit.tgc ? hit.tgc->phiCluster.hitList.front()->identify() : hit.prd->identify();
        
          if( m_doTruth ){
            if( !m_truthSummaryTool.empty() ) m_truthSummaryTool->add(id,1);
            if( truthHits.count(id) )       foundTruthHits.insert(id);
          }
          
          int nhits = hit.tgc ? hit.tgc->phiCluster.hitList.size() : 1;
          ATH_MSG_VERBOSE("findMaxima(Phi) phiHit " << m_idHelperSvc->toString(id) << " hits " << nhits );
        }
        
        maximum.sector = sector; // very fragile passing on of sector
        
        //check if the maximum is already filled, if so, don't add it again
        bool maximum_matched = false;
        for( auto pit = maxima.begin();pit!=maxima.end();++pit ){
          // reference to phi maximum
          MuonHough::MuonPhiLayerHough::Maximum& pmaximum = **pit;
          if (pmaximum.sector == maximum.sector && pmaximum.max == maximum.max && pmaximum.pos == maximum.pos && 
            pmaximum.hits.size() == maximum.hits.size() && pmaximum.binpos == maximum.binpos && 
            pmaximum.binposmin == maximum.binposmin && pmaximum.binposmax == maximum.binposmax){
            ATH_MSG_DEBUG("extendSeed: sector has already been added! Skip. ");
            bool maximum_hitmatched = true;//  check if there is a hit that is not the same
            for ( unsigned int k=0; k < maximum.hits.size(); ++k){
              if (maximum.hits[k] != pmaximum.hits[k]){// directly compare pointer address
                maximum_hitmatched = false;
                break;
              }
            }
            if (maximum_hitmatched){
              maximum_matched = true;
              break;
            }
          }
        }
        //remove the hits from hough
        hough.fillLayer2(maximum.hits,true);
        if (maximum_matched){
          //++nmaxima;
          continue;
        }
        else{
          maxima.push_back( new MuonHough::MuonPhiLayerHough::Maximum(maximum) );
          ++nmaxima;
        }
      }
      else{
        if( nmaxima > 0 ) {
          ATH_MSG_VERBOSE("findMaxima(Phi): No more maxima found " << nmaxima );
        }
        // ?!? same here, the function should return false if nothing was found, right?
        break;     
      }
    }
    hough.reset();
    return true;
  }

  void MuonLayerHoughTool::fillHitsPerSector(  std::set<Identifier>& truthHits,
                                               std::vector<std::unique_ptr<TgcHitClusteringObj>>& tgcClusteringObjs,
                                               const MuonLayerHoughTool::CollectionsPerSector& collectionsPerSector,
                                               const MdtPrepDataContainer*  mdtCont,  
                                               const CscPrepDataContainer*  /*cscCont*/,  
                                               const TgcPrepDataContainer*  tgcCont,  
                                               const RpcPrepDataContainer*  rpcCont,
                                               const sTgcPrepDataContainer* stgcCont,  
                                               const MMPrepDataContainer*   mmCont,
                                               MuonLayerHoughTool::HoughDataPerSector& houghData ) const {
    // loop over all possible station layers in the sector
    for( unsigned int tech=0;tech<m_ntechnologies;++tech ){
      for( unsigned int layerHash=0;layerHash<MuonStationIndex::sectorLayerHashMax();++layerHash ){
        const HashVec& hashes = collectionsPerSector.technologyRegionHashVecs[tech][layerHash];
        if( hashes.empty() ) continue;
        auto regionLayer = MuonStationIndex::decomposeSectorLayerHash( layerHash );
        
        HashVec::const_iterator iit = hashes.begin();
        HashVec::const_iterator iit_end = hashes.end();
        for( ;iit!=iit_end;++iit ){
          // !?! else if made by Felix
          if( mdtCont && mdtCont->size()>0 && tech == MuonStationIndex::MDT ) {
            auto pos = mdtCont->indexFindPtr(*iit);
            if( pos != nullptr ) fill(truthHits,*pos,houghData.hitVec[layerHash]);
          }
          else if( rpcCont && rpcCont->size()>0 && tech == MuonStationIndex::RPC ) {
            auto pos = rpcCont->indexFindPtr(*iit);
            if( pos != nullptr ) fill(truthHits,*pos,houghData.hitVec[layerHash],houghData.phiHitVec[regionLayer.first]);
          }
          else if( tgcCont && tgcCont->size()>0 && tech == MuonStationIndex::TGC ) {
            auto pos = tgcCont->indexFindPtr(*iit);
            if( pos != nullptr ) fill(truthHits, tgcClusteringObjs, *pos,houghData.hitVec[layerHash],houghData.phiHitVec[regionLayer.first],collectionsPerSector.sector);
          }
          else if( stgcCont && stgcCont->size()>0 && tech == MuonStationIndex::STGC ) {
            auto pos = stgcCont->indexFindPtr(*iit);
            if( pos != nullptr ) fill(truthHits,*pos,houghData.hitVec[layerHash],houghData.phiHitVec[regionLayer.first],collectionsPerSector.sector);
          }
          else if( mmCont && mmCont->size()>0 && tech == MuonStationIndex::MM ) {
            auto pos = mmCont->indexFindPtr(*iit);
            if( pos != nullptr ) fill(truthHits,*pos,houghData.hitVec[layerHash]);
          }
        }
      }
    }
  }

  void MuonLayerHoughTool::updateHits( MuonLayerHoughTool::PhiHitVec& hits, MuonHough::MuonPhiLayerHough& hough ) const {
    PhiHitVec::iterator hit = hits.begin();
    PhiHitVec::iterator hit_end = hits.end();
    for( ;hit!=hit_end;++hit ){
      if( (*hit)->debugInfo() ){
        float max = hough.maximum( (*hit)->r, (*hit)->phimin, (*hit)->phimax, (*hit)->debugInfo()->binpos );
        if( max > 100 ) ATH_MSG_WARNING(" Maximum value too large" << max );
        (*hit)->debugInfo()->ph = max;
        (*hit)->debugInfo()->rot = -99999.;
      }
      else{
        ATH_MSG_DEBUG("Failed to update hit: " << (*hit)->r << " " << (*hit)->phimin << " lay " << (*hit)->layer << " no debugInfo ");
      }
    }
  }

  void MuonLayerHoughTool::updateHits( MuonLayerHoughTool::HitVec& hits, MuonHough::MuonLayerHough& hough ) const {
    HitVec::iterator hit = hits.begin();
    HitVec::iterator hit_end = hits.end();
    for( ;hit!=hit_end;++hit ){
      if( (*hit)->debugInfo() ){
        std::pair<float,float> max = hough.maximum( (*hit)->x, (*hit)->ymin, (*hit)->debugInfo()->binpos, (*hit)->debugInfo()->bintheta );
        (*hit)->debugInfo()->ph = max.first;
        (*hit)->debugInfo()->rot = max.second;
      }
      else{
        ATH_MSG_DEBUG("Failed to update hit: " << (*hit)->x << " " << (*hit)->ymin << " lay " << (*hit)->layer << " no debugInfo ");
      }
    }
  }

  void MuonLayerHoughTool::matchTruth(std::set<Identifier>& truthHits, const PRD_MultiTruthCollection& truthCol, const Identifier& id, MuonHough::HitDebugInfo& debug ) const {
    typedef PRD_MultiTruthCollection::const_iterator iprdt;
    std::pair<iprdt, iprdt> range = truthCol.equal_range(id);
    // Loop over particles contributing to this cluster
    for(iprdt i = range.first; i != range.second; i++) {
      if(!i->second.isValid()) {
        ATH_MSG_WARNING("Unexpected invalid HepMcParticleLink in PRD_MultiTruthCollection");
      } else {
        const HepMcParticleLink& link = i->second;
        if( link.cptr() && abs(link.cptr()->pdg_id()) == 13 ){
          debug.barcode = link.barcode();
          debug.pdgId = link.cptr()->pdg_id();
          truthHits.insert(id);
        }
      }
    }
  }

  void MuonLayerHoughTool::fill( std::set<Identifier>& truthHits, const MdtPrepDataCollection& mdts, 
      MuonLayerHoughTool::HitVec& hits ) const {
    
    if( mdts.empty() ) return;
    auto truthCollections = m_truthNames.makeHandles();
    Identifier chid = mdts.identify();
    MuonStationIndex::DetectorRegionIndex region = m_idHelperSvc->regionIndex(chid);
    MuonStationIndex::LayerIndex layer = m_idHelperSvc->layerIndex(chid);
    int sector     = m_idHelperSvc->sector(chid);
    unsigned int technology = m_idHelperSvc->technologyIndex(chid);
    bool barrelLike = (region == MuonStationIndex::Barrel || layer == MuonStationIndex::BarrelExtended);
    unsigned int nmdts = 0;
    unsigned int nmdtsBad = 0;
    MdtPrepDataCollection::const_iterator mit = mdts.begin();
    MdtPrepDataCollection::const_iterator mit_end = mdts.end();
    for( ;mit!=mit_end;++mit ){
      const MdtPrepData& prd = **mit;
      if( prd.adc() < 50 || prd.status() != Muon::MdtStatusDriftTime ) {
        ++nmdtsBad;
        continue;
      }
      ++nmdts;
      const Identifier& id = prd.identify();

      float r = rCor(prd);
      float x = barrelLike ? r : prd.globalPosition().z();
      float y = barrelLike ? prd.globalPosition().z() : r;
      int sublayer = sublay(id);

      float ymin = y - prd.localPosition()[Trk::locR];
      float ymax = y + prd.localPosition()[Trk::locR];
      MuonHough::HitDebugInfo* debug = new MuonHough::HitDebugInfo(technology,sector,region,layer,sublayer);
      debug->time = prd.tdc();
      debug->r = prd.localPosition()[Trk::locR];
      
      std::map<unsigned int,unsigned int>::const_iterator pos = m_techToTruthNameIdx.find(technology);
      if (pos != m_techToTruthNameIdx.end()) {
        matchTruth(truthHits, *truthCollections[pos->second],id,*debug);
      }
      MuonHough::Hit* hit = new MuonHough::Hit(sublayer,x,ymin,ymax,1.,debug,&prd);
      hits.push_back(hit);
    }
    
    ATH_MSG_DEBUG("fillMDT: Filling " << m_idHelperSvc->toStringChamber(chid) << ": loc s"         << sector << " " << MuonStationIndex::regionName(region) << " " << MuonStationIndex::layerName(layer) 
                                                                           << " -> hits: "        << nmdts 
                                                                           << " bad "            << nmdtsBad 
                                                                           << " isSmallChamber " << m_idHelperSvc->isSmallChamber(chid) );


  }

  void MuonLayerHoughTool::fill( std::set<Identifier>& truthHits, const RpcPrepDataCollection& rpcs, 
      MuonLayerHoughTool::HitVec& hits, MuonLayerHoughTool::PhiHitVec& phiHits ) const {
    
    if( rpcs.empty() ) return;
    auto truthCollections = m_truthNames.makeHandles();
    Identifier chid = rpcs.identify();
    unsigned int technology                      = m_idHelperSvc->technologyIndex(chid);
    MuonStationIndex::LayerIndex layer           = m_idHelperSvc->layerIndex(chid);
    MuonStationIndex::DetectorRegionIndex region = m_idHelperSvc->regionIndex(chid);
    int sector                                   = m_idHelperSvc->sector(chid);

    // check whether there are eta and phi hits
    unsigned int neta = 0;
    unsigned int nphi = 0;
    RpcPrepDataCollection::const_iterator mit = rpcs.begin();
    RpcPrepDataCollection::const_iterator mit_end = rpcs.end();
    for( ;mit!=mit_end;++mit ){
      if( m_idHelperSvc->rpcIdHelper().measuresPhi((*mit)->identify()) )  ++nphi;
      else ++neta;
    }
    ATH_MSG_DEBUG("fillRPC: Filling " << m_idHelperSvc->toStringChamber(chid) 
                                      << ": loc s"    << sector << " " << MuonStationIndex::regionName(region) << " " << MuonStationIndex::layerName(layer) 
                                      << " -> eta hits " << neta 
                                      << " phi hits "  << nphi );


    mit = rpcs.begin();
    for( ;mit!=mit_end;++mit ){
      const RpcPrepData& prd = **mit;
      const Identifier& id = prd.identify();
      int sublayer = sublay(id);
      MuonHough::HitDebugInfo* debug = new MuonHough::HitDebugInfo(technology,sector,region,layer,sublayer);
      debug->isEtaPhi = (neta && nphi);
      debug->trigConfirm = 1;
      debug->time = prd.time();
      std::map<unsigned int,unsigned int>::const_iterator pos = m_techToTruthNameIdx.find(technology);
      if (pos != m_techToTruthNameIdx.end()) {
        matchTruth(truthHits, *truthCollections[pos->second],id,*debug);
      }
      float weight = (neta && nphi) ? 2 : 1;
      if( m_idHelperSvc->rpcIdHelper().measuresPhi(id) ) {
        const float r = rCor(prd);
        const float phi = prd.globalPosition().phi();
        const double phi1 = phi; //phiCor(phi,sector);
        debug->r = -99999;
        MuonHough::PhiHit* hit = new MuonHough::PhiHit(sublayer,r,phi1,phi1,weight,debug,&prd);
        phiHits.push_back(hit);
      }
      else{
        const float x = rCor(prd);
        const float y = prd.globalPosition().z();
        const float stripCor = 0.5*prd.detectorElement()->StripWidth(false);
        const float ymin = y - stripCor;
        const float ymax = y + stripCor;
        debug->r = stripCor;
        MuonHough::Hit* hit = new MuonHough::Hit(sublayer,x,ymin,ymax,weight,debug,&prd);
        hits.push_back(hit);
      }
    }
  }

  void MuonLayerHoughTool::fill( std::set<Identifier>& truthHits, const MMPrepDataCollection& mms, 
      MuonLayerHoughTool::HitVec& hits ) const {
    
    if( mms.empty() ) return;
    auto truthCollections = m_truthNames.makeHandles();
    Identifier chid = mms.identify();
    MuonStationIndex::DetectorRegionIndex region = m_idHelperSvc->regionIndex(chid);
    MuonStationIndex::LayerIndex layer = m_idHelperSvc->layerIndex(chid);
    int sector = m_idHelperSvc->sector(chid);
    unsigned int technology = m_idHelperSvc->technologyIndex(chid);
    ATH_MSG_DEBUG("fillMM: Filling " << m_idHelperSvc->toStringChamber(chid) << ": loc s" << sector << " " << MuonStationIndex::regionName(region) << " " << MuonStationIndex::layerName(layer) 
                                                                          << " -> hits "  << mms.size() );

    MMPrepDataCollection::const_iterator mit = mms.begin();
    MMPrepDataCollection::const_iterator mit_end = mms.end();
    for( ;mit!=mit_end;++mit ){
      const MMPrepData& prd = **mit;
      const Identifier& id = prd.identify();

      float x = prd.globalPosition().z();
      float y = rCor(prd);
      int sublayer = sublay(id);
      float stripCor = 0.45; // get from det el
      float ymin = y - stripCor;
      float ymax = y + stripCor;
      MuonHough::HitDebugInfo* debug = new MuonHough::HitDebugInfo(technology,sector,region,layer,sublayer);
      debug->r = stripCor;
      std::map<unsigned int,unsigned int>::const_iterator pos = m_techToTruthNameIdx.find(technology);
      if (pos != m_techToTruthNameIdx.end()) {
        matchTruth(truthHits, *truthCollections[pos->second],id,*debug);
      }
      MuonHough::Hit* hit = new MuonHough::Hit(sublayer,x,ymin,ymax,1.,debug,&prd);
      hits.push_back(hit);
    }
  }

  void MuonLayerHoughTool::fill( std::set<Identifier>& truthHits, const sTgcPrepDataCollection& stgcs, 
      MuonLayerHoughTool::HitVec& hits, MuonLayerHoughTool::PhiHitVec& phiHits, int selectedSector ) const {
    
    if( stgcs.empty() ) return;
    auto truthCollections = m_truthNames.makeHandles();
    Identifier chid = stgcs.identify();
    MuonStationIndex::DetectorRegionIndex region = m_idHelperSvc->regionIndex(chid);
    MuonStationIndex::LayerIndex layer = m_idHelperSvc->layerIndex(chid);
    int sector = m_idHelperSvc->sector(chid);
    bool isNeighbouringSector = sector!=selectedSector;
    unsigned int technology = m_idHelperSvc->technologyIndex(chid);
    ATH_MSG_DEBUG("fillsTGC: Filling " << m_idHelperSvc->toStringChamber(chid) 
                         << ": loc s" << sector << " " << MuonStationIndex::regionName(region) << " " << MuonStationIndex::layerName(layer) 
                         << " -> hits: " << stgcs.size() );

    sTgcPrepDataCollection::const_iterator mit = stgcs.begin();
    sTgcPrepDataCollection::const_iterator mit_end = stgcs.end();
    for( ;mit!=mit_end;++mit ){
      const sTgcPrepData& prd = **mit;
      const Identifier& id = prd.identify();
      int channelType = m_idHelperSvc->stgcIdHelper().channelType(id);

      // only pick up phi hits in neighbouring sectors
      if( isNeighbouringSector && channelType == 1 ) continue;
      if( m_onlyUseCurrentBunch && (prd.getBcBitMap() & sTgcPrepData::BCBIT_CURRENT) != sTgcPrepData::BCBIT_CURRENT ) continue;
      int sublayer = sublay(id);
      
      MuonHough::HitDebugInfo* debug = new MuonHough::HitDebugInfo(technology,sector,region,layer,sublayer);
      debug->isEtaPhi = 1;
      debug->trigConfirm = (prd.getBcBitMap() & sTgcPrepData::BCBIT_CURRENT) == sTgcPrepData::BCBIT_CURRENT;
      debug->time = prd.getBcBitMap();

      std::map<unsigned int,unsigned int>::const_iterator pos = m_techToTruthNameIdx.find(technology);
      if (pos != m_techToTruthNameIdx.end()) {
        matchTruth(truthHits, *truthCollections[pos->second],id,*debug);
      }
      if( m_idHelperSvc->stgcIdHelper().channelType(id) == 1 ) {
        // eta strips
        float x = prd.globalPosition().z();
        float y = rCor(prd);
        float stripCor = 1.5; // get from det el
        const MuonGM::MuonChannelDesign* design = prd.detectorElement()->getDesign( id );
        if(design) {
          double stripWidth=design->inputWidth;
          double stripLength=design->channelLength(m_idHelperSvc->stgcIdHelper().channel(id));
          if(m_debugHough) ATH_MSG_DEBUG(" eta strip width " << stripWidth << " stripLength " << stripLength);
          stripCor = 0.5*stripWidth;
        }
        debug->r = stripCor;
        float ymin = y - stripCor;
        float ymax = y + stripCor;
        MuonHough::Hit* hit = new MuonHough::Hit(sublayer,x,ymin,ymax,1.,debug,&prd);
        hits.push_back(hit);
      }
      else{
        double chWidth = 0;
        if( m_idHelperSvc->stgcIdHelper().channelType(id) == 0 ) {
            
          // pads
          const MuonGM::MuonPadDesign* design = prd.detectorElement()->getPadDesign(id);
          if( !design ) {
            ATH_MSG_WARNING("No design found for " << m_idHelperSvc->toString(id) );
            delete debug;
            continue;
          }
          // weird large number
          double chWidthOLD = 0.5*design->channelWidth(prd.localPosition(),true);

          // inputPhiPitch is in degrees
          float radius = prd.globalPosition().perp();
          chWidth = 0.5*design->inputPhiPitch*M_PI/180.*radius;
          if(m_debugHough) ATH_MSG_DEBUG(" sPadWidth " << design->sPadWidth  << " lPadWidth " << design->lPadWidth  << " inputRowWidth " << design->inputRowWidth);

          if(m_debugHough) ATH_MSG_DEBUG(" Pad chWidth " << chWidth  << " OLD " << chWidthOLD << " phi global " << prd.globalPosition().phi());
        }
        else if( m_idHelperSvc->stgcIdHelper().channelType(id) == 2 ) {
          const MuonGM::MuonChannelDesign* design = prd.detectorElement()->getDesign(id);
          if( !design ) {
            ATH_MSG_WARNING("No design found for " << m_idHelperSvc->toString(id) );
            delete debug;
            continue;
          }
          // double etaWidth=design->channelLength(idhelper->channel(id));
          double phiMaxWidth=design->maxYSize/design->nch;
          // double phiMinWidth=design->minYSize/design->nch;
          double chWidthOLD = 0.5*design->channelWidth(prd.localPosition());
          chWidth = 0.5*phiMaxWidth;

          if(m_debugHough) ATH_MSG_DEBUG(" Wire Gang chWidth " << chWidth  << " OLD " << chWidthOLD << " phi global " << prd.globalPosition().phi());
        }
        
        Amg::Vector2D lp1(prd.localPosition().x()+chWidth,prd.localPosition().y());
        Amg::Vector3D gp1;
        prd.detectorElement()->surface(id).localToGlobal(lp1,gp1,gp1);
        
        lp1[0] = prd.localPosition().x()-chWidth;
        Amg::Vector3D gp2;
        prd.detectorElement()->surface(id).localToGlobal(lp1,gp2,gp2);
        
        double phi1 = gp1.phi();
        double phi2 = gp2.phi();
        double phi1c = phi1;//phiCor(phi1,selectedSector);
        double phi2c = phi2;//phiCor(phi2,selectedSector);
        if( fabs(phi1c-phi2c) > 0.3 ){
          ATH_MSG_WARNING("bad local phi: in " << phi1 << ", " << phi2 << " sector phi " << m_sectorMapping.sectorPhi(selectedSector) << " phicor " << phi1c << ", " << phi2c );
        }
        if( isNeighbouringSector && !(m_sectorMapping.insideSector(selectedSector,phi1)||m_sectorMapping.insideSector(selectedSector,phi2)) ){
          ATH_MSG_DEBUG("Dropping phi hit in neighbouring sector " << m_idHelperSvc->toString(id) << " phi min " 
                        << std::min(phi1c,phi2c) << " max " << std::max(phi1c,phi2c)
                        << " global phi: in " << phi1 << ", " << phi2 << " sector phi " << m_sectorMapping.sectorPhi(selectedSector) );
          delete debug;
          continue;
        } 
        float r = rCor(prd);
        MuonHough::PhiHit* phiHit = new MuonHough::PhiHit(sublayer,r,std::min(phi1c,phi2c),std::max(phi1c,phi2c),1,debug,&prd);
        ATH_MSG_VERBOSE("Phi hit " << m_idHelperSvc->toString(id) << " r " << r << " phi min " << phiHit->phimin << " phi max " << phiHit->phimax 
            << " bc " << debug->barcode << " chw " << chWidth << " trigC " << debug->trigConfirm << " g phi " << phi1 << " " << phi2 );
        phiHits.push_back(phiHit);
      }
    }
  }

  void MuonLayerHoughTool::fill( std::set<Identifier>& truthHits, std::vector<std::unique_ptr<TgcHitClusteringObj>>& tgcClusteringObjs, 
      const TgcPrepDataCollection& tgcs, MuonLayerHoughTool::HitVec& hits, MuonLayerHoughTool::PhiHitVec& phiHits, 
      int sector ) const {
    
    if( tgcs.empty() ) return;
    tgcClusteringObjs.push_back( std::make_unique<TgcHitClusteringObj>(&m_idHelperSvc->tgcIdHelper()) );
    TgcHitClusteringObj& clustering = *tgcClusteringObjs.back();
    std::vector<const TgcPrepData*> prds;
    prds.insert(prds.begin(),tgcs.begin(),tgcs.end());
    clustering.cluster(prds);
    clustering.buildClusters3D();

    Identifier chid = tgcs.identify();
    MuonStationIndex::DetectorRegionIndex region = m_idHelperSvc->regionIndex(chid);
    MuonStationIndex::LayerIndex layer = m_idHelperSvc->layerIndex(chid);

    if( clustering.clusters3D.empty() ) {
      ATH_MSG_DEBUG("TgcHitClusteringObj, no 3D clusters! ");
      if( msgLvl(MSG::DEBUG) ){
        for(std::vector<const TgcPrepData*>::iterator it=prds.begin();it!=prds.end();++it ){
          ATH_MSG_DEBUG("   " << m_idHelperSvc->toString( (*it)->identify() ));
        }
      }
      return;
    }
    if( !clustering.bestEtaCluster() || clustering.bestEtaCluster()->hitList.empty() || !clustering.bestEtaCluster()->hitList.front()){
      ATH_MSG_DEBUG("TgcHitClusteringObj, no eta cluster selected! ");
      if( msgLvl(MSG::DEBUG) ){
        for(std::vector<const TgcPrepData*>::iterator it=prds.begin();it!=prds.end();++it ){
          ATH_MSG_DEBUG("   " << m_idHelperSvc->toString( (*it)->identify() ));
        }
      }
      return;
    }
    auto truthCollections = m_truthNames.makeHandles();
    std::vector<int> sectors;
    getSectors(clustering.clusters3D.front(),sectors);
    unsigned int technology = m_idHelperSvc->technologyIndex(chid);
    for( unsigned int si=0;si<sectors.size();++si ){
      if( sectors[si] != sector ) continue;
      std::vector< TgcClusterObj3D >::const_iterator cl_it = clustering.clusters3D.begin();
      std::vector< TgcClusterObj3D >::const_iterator cl_it_end = clustering.clusters3D.end();
      
      for( ;cl_it!=cl_it_end;++cl_it ){
  const TgcClusterObj3D& cl = *cl_it;
  if( cl.etaCluster.hitList.empty() ) {
    ATH_MSG_WARNING("Incomplete TgcClusterObj3D in chamber " << m_idHelperSvc->toString(chid) );
    continue;
  }
  const Identifier& id = cl.etaCluster.hitList.front()->identify();

  double x = cl.p11.z();
  double y11 = rCor( cl, 1, sector );
  double y12 = rCor( cl, 2, sector );
  double y21 = rCor( cl, 3, sector );
  double y22 = rCor( cl, 4, sector );
  double phi11 = cl.p11.phi();
  double phi12 = cl.p12.phi();
  double phi21 = cl.p21.phi();
  double phi22 = cl.p22.phi();
  double ymin = std::min( std::min(y11,y12), std::min(y21,y22) );
  double ymax = std::max( std::max(y11,y12), std::max(y21,y22) );
  double phimin = std::min( std::min(phi11,phi12), std::min(phi21,phi22) );
  double phimax = std::max( std::max(phi11,phi12), std::max(phi21,phi22) );
  double phi1 = phimin; //phiCor(phimin,sector);
  double phi2 = phimax; //phiCor(phimax,sector);
  int sublayer = sublay(id,x);

        MuonHough::HitDebugInfo* debug = new MuonHough::HitDebugInfo(technology,sector,region,layer,sublayer);
  debug->clusterSize = cl.etaCluster.hitList.size();
  debug->clusterLayers = cl.etaCluster.layers();
  debug->isEtaPhi = cl.phiCluster.layers();
  debug->time = cl.etaCluster.hitList.front()->getBcBitMap();
  std::map<unsigned int,unsigned int>::const_iterator pos = m_techToTruthNameIdx.find(technology);
  if (pos != m_techToTruthNameIdx.end()) {
    matchTruth(truthHits, *truthCollections[pos->second],id,*debug);
  }

  MuonHough::HitDebugInfo* phiDebug = new MuonHough::HitDebugInfo(*debug);
  phiDebug->clusterSize = cl.phiCluster.hitList.size();
  phiDebug->clusterLayers = cl.phiCluster.layers();
  phiDebug->isEtaPhi = cl.etaCluster.layers();

  MuonHough::Hit* hit = new MuonHough::Hit(sublayer,x,ymin,ymax,2*cl.etaCluster.layers(),debug,0,&cl);
  MuonHough::PhiHit* phiHit = new MuonHough::PhiHit(sublayer,y11,phi1,phi2,2*cl.phiCluster.layers(), phiDebug,0,&cl);
  hits.push_back(hit);
  phiHits.push_back(phiHit);
      }
    }
    ATH_MSG_DEBUG("fillTGC: Filling " << m_idHelperSvc->toStringChamber(chid) << ": loc s"    << sector << " " << MuonStationIndex::regionName(region) << " " << MuonStationIndex::layerName(layer) 
                                                                           << " -> etaHits: " << hits.size() 
                                                                           << " phiHits: "  << phiHits.size() 
                                                                           << " sectors: "  << sectors.size() );
  }
 
  void MuonLayerHoughTool::insertHash( const IdentifierHash& hash, const Identifier& id ) {
    insertHash(m_idHelperSvc->sector(id),hash,id);
  }

  void MuonLayerHoughTool::insertHash( int sector, const IdentifierHash& hash, const Identifier& id ) {
    MuonStationIndex::TechnologyIndex techIndex = m_idHelperSvc->technologyIndex(id);
    int sectorLayerHash = MuonStationIndex::sectorLayerHash(m_idHelperSvc->regionIndex(id),m_idHelperSvc->layerIndex(id));
    m_collectionsPerSector[sector-1].technologyRegionHashVecs[techIndex][sectorLayerHash].push_back(hash);
  }

  // all chambers are mapped onto a layer and sector map
  void MuonLayerHoughTool::initializeSectorMapping(const MuonGM::MuonDetectorManager* detMgr) {
    m_collectionsPerSector.resize(MuonStationIndex::numberOfSectors());
    // set sector numbers
    unsigned int nsectorHashMax = MuonStationIndex::sectorLayerHashMax();
    for( unsigned int i=0;i<m_collectionsPerSector.size();++i ) {
      m_collectionsPerSector[i].sector=i+1;
      m_collectionsPerSector[i].technologyRegionHashVecs.resize(m_ntechnologies);
      for( auto it = m_collectionsPerSector[i].technologyRegionHashVecs.begin();it!=m_collectionsPerSector[i].technologyRegionHashVecs.end(); ++it ) {
        it->resize(nsectorHashMax);
      }
    }
    ATH_MSG_DEBUG("Initializing hashes: number of sectors " << MuonStationIndex::numberOfSectors() 
                  << " technologies " << m_ntechnologies << " sectorLayers " << MuonStationIndex::sectorLayerHashMax() );
    // loop over all available MDT collection identifiers and order them per sector
    MuonIdHelper::const_id_iterator it = m_idHelperSvc->mdtIdHelper().module_begin();
    MuonIdHelper::const_id_iterator it_end = m_idHelperSvc->mdtIdHelper().module_end();
    for( ;it!=it_end; ++it ){
      IdentifierHash hash;
      m_idHelperSvc->mdtIdHelper().get_module_hash(*it,hash);
      insertHash(hash,*it);
    }

    // loop over all available RPC collection identifiers and order them per sector
    it = m_idHelperSvc->rpcIdHelper().module_begin();
    it_end = m_idHelperSvc->rpcIdHelper().module_end();
    for( ;it!=it_end; ++it ){
      IdentifierHash hash;
      m_idHelperSvc->rpcIdHelper().get_module_hash(*it,hash);
      insertHash(hash,*it);
    }

    // loop over all available CSC collection identifiers and order them per sector
    if (m_idHelperSvc->hasCSC()) {
      it = m_idHelperSvc->cscIdHelper().module_begin();
      it_end = m_idHelperSvc->cscIdHelper().module_end();
      for( ;it!=it_end; ++it ){
        IdentifierHash hash;
        m_idHelperSvc->cscIdHelper().get_module_hash(*it,hash);
        insertHash(hash,*it);
      }
    }

    // loop over all available MM collection identifiers and order them per sector
    if (m_idHelperSvc->hasMM()) {
      it = m_idHelperSvc->mmIdHelper().detectorElement_begin();
      it_end = m_idHelperSvc->mmIdHelper().detectorElement_end();
      for( ;it!=it_end; ++it ){
        IdentifierHash hash;
        m_idHelperSvc->mmIdHelper().get_module_hash(*it,hash);
        insertHash(hash,*it);
      }
    }

    // loop over all available STGC collection identifiers and order them per sector
    if (m_idHelperSvc->hasSTgc()) {
      it = m_idHelperSvc->stgcIdHelper().detectorElement_begin();
      it_end = m_idHelperSvc->stgcIdHelper().detectorElement_end();
      for( ;it!=it_end; ++it ){
        IdentifierHash hash;
        m_idHelperSvc->stgcIdHelper().get_module_hash(*it,hash);
        int sector = m_idHelperSvc->sector(*it);
        insertHash(sector,hash,*it);
        int sectorU = sector != 1 ? sector-1 : 16;
        int sectorD = sector != 16 ? sector+1 : 1;
        insertHash(sectorU,hash,*it);
        insertHash(sectorD,hash,*it);
      }
    }

    // loop over all available TGC collection identifiers and order them per sector
    it = m_idHelperSvc->tgcIdHelper().module_begin();
    it_end = m_idHelperSvc->tgcIdHelper().module_end();
    for( ;it!=it_end; ++it ){
      const MuonGM::TgcReadoutElement* detEl = detMgr->getTgcReadoutElement(*it);
      if( !detEl ) {
        ATH_MSG_DEBUG(" No detector element found for " << m_idHelperSvc->toString(*it) );
        continue;
      }
      IdentifierHash hash;
      m_idHelperSvc->tgcIdHelper().get_module_hash(*it,hash);
      int nstrips = detEl->getNStrips(1);
      Amg::Vector3D p1 = detEl->channelPos(1,1,1);
      Amg::Vector3D p2 = detEl->channelPos(1,1,nstrips);
      std::vector<int> sectors1;
      getSectors(p1,sectors1);
      std::set<int> added;
      std::vector<int>::iterator sit = sectors1.begin();
      std::vector<int>::iterator sit_end = sectors1.end();
      for( ;sit!=sit_end; ++sit ){
        insertHash(*sit,hash,*it);
        added.insert(*sit);
      }

      std::vector<int> sectors2;
      getSectors(p2,sectors2);
      sit = sectors2.begin();
      sit_end = sectors2.end();
      for( ;sit!=sit_end; ++sit ){
        if( added.count(*sit) ) continue;
        added.insert(*sit);
        insertHash(*sit,hash,*it);
      }

    }

    if( msgLvl(MSG::DEBUG) ) ATH_MSG_DEBUG(" Printing collections per sector, number of technologies " << m_ntechnologies);
    for( int sector = 1; sector<=16; ++sector ){
      MuonStationIndex::DetectorRegionIndex currentRegion = MuonStationIndex::DetectorRegionUnknown;
      if( msgLvl(MSG::DEBUG) ) ATH_MSG_DEBUG(" sector " << sector);
      TechnologyRegionHashVec& vec = m_collectionsPerSector[sector-1].technologyRegionHashVecs;
      for( unsigned int hash = 0; hash < nsectorHashMax; ++hash ){
        std::pair<MuonStationIndex::DetectorRegionIndex,MuonStationIndex::LayerIndex> regionLayer = MuonStationIndex::decomposeSectorLayerHash(hash);
        if( msgLvl(MSG::DEBUG) ) if( regionLayer.first != currentRegion ) ATH_MSG_DEBUG("  " << MuonStationIndex::regionName(regionLayer.first));
        bool first = true;
        currentRegion = regionLayer.first;
        for( unsigned int tech=0; tech<m_ntechnologies;++tech ){
          std::stable_sort(vec[tech][hash].begin(),vec[tech][hash].end());
          if( !vec[tech][hash].empty() ) {
            if( msgLvl(MSG::DEBUG) ) {
              if( first ) {
                ATH_MSG_DEBUG("  " << std::setw(7) << MuonStationIndex::layerName(regionLayer.second));
                first = false;
              }
              ATH_MSG_DEBUG(" " << std::setw(4) << MuonStationIndex::technologyName(static_cast<MuonStationIndex::TechnologyIndex>(tech)) 
                              << " " << std::setw(4) << vec[tech][hash].size()); 
            }
          }
        }
      }
    }
  }

  void MuonLayerHoughTool::printTruthSummary( std::set<Identifier>& truth, std::set<Identifier>& found ) const {
    if( truth.size() == found.size() ){
      ATH_MSG_DEBUG(" All hits found: truth " << truth.size() << " found " << found.size() );
    }
    else{
      ATH_MSG_DEBUG(" Some truth hits not found: truth " << truth.size() << " found " << found.size() );
      std::vector<Identifier> result(truth.size()-found.size());
      std::vector<Identifier>::iterator pos = std::set_difference(truth.begin(),truth.end(),found.begin(),found.end(),result.begin());
      result.resize(pos-result.begin());
      for( std::vector<Identifier>::iterator it=result.begin();it!=result.end();++it ){
        ATH_MSG_DEBUG("  " << m_idHelperSvc->toString(*it) );
      }
    }
  }

  void MuonLayerHoughTool::fillNtuple( MuonLayerHoughTool::HoughDataPerSectorVec& houghDataPerSectorVec ) const {
    for( const auto& it : houghDataPerSectorVec.vec ){
      for( const auto& rit : it.hitVec ) {
        m_ntuple->fill(rit);
      }
      for( const auto& rit : it.phiHitVec ) {
        m_ntuple->fill(rit);
      } 
    }
  }
}
