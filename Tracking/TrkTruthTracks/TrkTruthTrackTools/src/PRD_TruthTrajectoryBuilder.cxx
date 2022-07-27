/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////
// PRD_TruthTrajectoryBuilder.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

// package include
#include "PRD_TruthTrajectoryBuilder.h"
// Trk
#include "TrkToolInterfaces/IPRD_Provider.h"
#include "TrkTruthTrackInterfaces/IPRD_TruthTrajectoryManipulator.h"
#include "TrkEventUtils/PrepRawDataComparisonFunction.h"
// DetectorDescription
#include "AtlasDetDescr/AtlasDetectorID.h"
// HepMC
#include "AtlasHepMC/GenParticle.h"
#include "AtlasHepMC/GenVertex.h"


/** Constructor **/
Trk::PRD_TruthTrajectoryBuilder::PRD_TruthTrajectoryBuilder(const std::string& t, const std::string& n, const IInterface* p) : 
  AthAlgTool(t,n,p),
  m_idHelper(nullptr)
{
    declareInterface<Trk::IPRD_TruthTrajectoryBuilder>(this);
    // the PRD providers that turn Identifier -> IdentiferHash and get the PRD
}

// Athena algtool's Hooks - initialize
StatusCode  Trk::PRD_TruthTrajectoryBuilder::initialize()
{
    std::cout << "Nora: min pT " << m_minPt << "\n";
    
    ATH_MSG_INFO("Initializing ...");
    // Set up ATLAS ID helper to be able to identify the PRD's det-subsystem
    if (detStore()->retrieve(m_idHelper, "AtlasID").isFailure()) {
         ATH_MSG_ERROR ("Could not get AtlasDetectorID helper. Arborting ...");
         return StatusCode::FAILURE;
    }
    // get the ID PRD Provider
    if ( !m_idPrdProvider.empty() && m_idPrdProvider.retrieve().isFailure()){
        ATH_MSG_ERROR ("Could not get " << m_idPrdProvider << ". Arborting ..." );
        return StatusCode::FAILURE;
    }
    // get the Muon System PRD Provider
    if ( !m_msPrdProvider.empty() && m_msPrdProvider.retrieve().isFailure()){
        ATH_MSG_ERROR ("Could not get " << m_msPrdProvider << ". Arborting ..." );
        return StatusCode::FAILURE;
    }     
    // get the manipulators
    if ( !m_prdTruthTrajectoryManipulators.empty()  && m_prdTruthTrajectoryManipulators.retrieve().isFailure()){
        ATH_MSG_ERROR ("Could not get configured " << m_prdTruthTrajectoryManipulators << ". Arborting ..." );
        return StatusCode::FAILURE;
    }
    
    ATH_CHECK( m_prdMultiTruthCollectionNames.initialize() );

    return StatusCode::SUCCESS;
}

StatusCode Trk::PRD_TruthTrajectoryBuilder::refreshEvent()  {

   ATH_MSG_INFO("Calling refreshEvent() to reset cache and retrieve collections");
   // clear the cache & reserve
   m_gpPrdTruthTrajectories.clear();
   m_prdMultiTruthCollections.clear();
   m_prdMultiTruthCollections.reserve(m_prdMultiTruthCollectionNames.size());
   // load the PRD collections from SG
   for(const auto& pmtCollNameIter:m_prdMultiTruthCollectionNames){
     // try to retrieve the PRD multi truth collection
     SG::ReadHandle<PRD_MultiTruthCollection> curColl (pmtCollNameIter);
     if (!curColl.isValid()){
       ATH_MSG_WARNING("Could not retrieve " << pmtCollNameIter << ". Ignoring ... ");
     }
     else{
       ATH_MSG_INFO("Added " << pmtCollNameIter << " to collection list for truth track creation.");
       m_prdMultiTruthCollections.push_back(curColl.cptr());
     }
   }
   // retrieve collection call to the PRD_Providers
   if (!m_idPrdProvider.empty() && m_idPrdProvider->retrieveCollection().isFailure()){
       ATH_MSG_ERROR ("Failure in collection retrieval of " << m_idPrdProvider << ". Arborting ..." );
       return StatusCode::FAILURE; 
   }
   if (!m_msPrdProvider.empty() && m_msPrdProvider->retrieveCollection().isFailure()){
       ATH_MSG_ERROR ("Failure in collection retrieval of " << m_msPrdProvider << ". Arborting ..." );
       return StatusCode::FAILURE; 
   }
   // all good
   return StatusCode::SUCCESS;
   
}

const std::map<HepMC::ConstGenParticlePtr, Trk::PRD_TruthTrajectory >& Trk::PRD_TruthTrajectoryBuilder::truthTrajectories() const {
    // ndof
    size_t ndofTotal = 0;
    size_t ndof      = 0;
    // PART 1 --------------------------------------------------------------------------------------------------------
    // loop over the PRD_MultiTruthCollection, search for the PRD and create (if necessary and entry in the return map)
    std::vector<const PRD_MultiTruthCollection*>::const_iterator pmtCollIter  = m_prdMultiTruthCollections.begin();
    std::vector<const PRD_MultiTruthCollection*>::const_iterator pmtCollIterE = m_prdMultiTruthCollections.end();
    for ( ; pmtCollIter != pmtCollIterE; ++pmtCollIter ){
        // loop over the map and get the identifier, GenParticle relation
        std::cout << "Nora: size of collection: " << (*pmtCollIter)->size() << "\n";
        PRD_MultiTruthCollection::const_iterator prdMtCIter  = (*pmtCollIter)->begin();
        PRD_MultiTruthCollection::const_iterator prdMtCIterE = (*pmtCollIter)->end();
        for ( ; prdMtCIter != prdMtCIterE; ++ prdMtCIter ){

            // check if entry exists and if   
#ifdef HEPMC3
            HepMC::ConstGenParticlePtr curGenP       = (*prdMtCIter).second.scptr();
#else
//AV Looks like an implicit conversion
            HepMC::ConstGenParticlePtr curGenP       = (*prdMtCIter).second;
#endif
            Identifier                curIdentifier = (*prdMtCIter).first;
            if(m_idHelper->is_trt(curIdentifier)) std::cout << "I'm a trt thingy hi\n";
            // apply the min pT cut 
            std::cout << "Nora: pt " << curGenP->momentum().perp() << " pgid " << std::abs(curGenP->pdg_id()) << " geantino " << m_geantinos << "\n";
            if ( curGenP->momentum().perp() < m_minPt ) continue;
            // skip geantinos if required
            if (!m_geantinos && std::abs(curGenP->pdg_id())==999) continue;
            // get the associated PRD from the provider
            std::cout << "Is id? " << m_idHelper->is_indet(curIdentifier) << "\n";
            const Trk::PrepRawData* prd = m_idHelper->is_indet(curIdentifier) ?
                m_idPrdProvider->prdFromIdentifier(curIdentifier,ndof) : m_msPrdProvider->prdFromIdentifier(curIdentifier,ndof);
            // stuff it into the trajectory if you found a PRD
            if (prd){
                // try to find the entry for this GenParticle 
                auto prdTrajIter = m_gpPrdTruthTrajectories.find(curGenP);
                if ( prdTrajIter ==  m_gpPrdTruthTrajectories.end() ){
                    // first PRD associated to this: create PRD_TruthTrajectory object
                    Trk::PRD_TruthTrajectory newPrdTruthTrajectory;
                    newPrdTruthTrajectory.prds.push_back(prd);
                    newPrdTruthTrajectory.nDoF = ndof-5;
                    // register the GenParticle only once
                    newPrdTruthTrajectory.genParticle = curGenP;
                    // fill into map
                    m_gpPrdTruthTrajectories[curGenP] = newPrdTruthTrajectory;
                    ndofTotal = ndof;
                } else {
                    // this PRD_TruthTrajectory already exists
                    (prdTrajIter->second).prds.push_back(prd);
                    (prdTrajIter->second).nDoF += ndof;
                    ndofTotal = (prdTrajIter->second).nDoF;
                }
                ATH_MSG_INFO("  Associating PRD with " << ndof << " degrees of freedom, total N.d.o.F : " << ndofTotal );
                ATH_MSG_INFO("  Associating Identifier " << curIdentifier << " with particle at [ " << curGenP << " ]." );
                std::string prdtype = m_idHelper->is_pixel(curIdentifier) ? "Pixel" : m_idHelper->is_sct(curIdentifier) ?  "SCT" : "TRT";
                ATH_MSG_INFO("  PRD is a " << prdtype);
            } else {
                std::string prdtype = m_idHelper->is_pixel(curIdentifier) ? "Pixel" : m_idHelper->is_sct(curIdentifier) ?  "SCT" : "TRT";
                ATH_MSG_INFO("  Failed to get " << prdtype << " PRD");
            }
        }        
    }
    std::cout << "Nora found " << m_gpPrdTruthTrajectories.size() << " tracks\n";
    // PART 2 --------------------------------------------------------------------------------------------------------
    // loop through the provided list of manipulators ( sorter is included )
    auto prdTruthTrajIter  = m_gpPrdTruthTrajectories.begin();
    auto prdTruthTrajIterE = m_gpPrdTruthTrajectories.end();
    for ( ; prdTruthTrajIter != prdTruthTrajIterE; ++prdTruthTrajIter ){
        if ( !m_prdTruthTrajectoryManipulators.empty() ){
            ToolHandleArray<IPRD_TruthTrajectoryManipulator>::const_iterator prdTTMIter  = m_prdTruthTrajectoryManipulators.begin();
            ToolHandleArray<IPRD_TruthTrajectoryManipulator>::const_iterator prdTTMIterE = m_prdTruthTrajectoryManipulators.end();
            for ( ; prdTTMIter != prdTTMIterE; ++prdTTMIter ){
                if ((*prdTTMIter)->manipulateTruthTrajectory((*prdTruthTrajIter).second))
                    ATH_MSG_INFO("PRD truth trajectory got manipulated by: " << (*prdTTMIter).name() );
            }
        }
    }
    // return the truth trajectories and leave it to the TruthTrack creation to proceed further
    return m_gpPrdTruthTrajectories;
}
                                    
StatusCode  Trk::PRD_TruthTrajectoryBuilder::finalize()
{
    // clear the cache a last time
    m_gpPrdTruthTrajectories.clear();
    m_prdMultiTruthCollections.clear();    
    ATH_MSG_INFO("Finalizing ...");
    return StatusCode::SUCCESS;
}

