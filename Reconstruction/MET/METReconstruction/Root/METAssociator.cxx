///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// METAssociator.cxx
// Implementation for class METAssociator
//
// This is the base class for tools that construct MET terms
// from other object collections.
//
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//
// Author: P Loch, S Resconi, TJ Khoo, AS Mete
///////////////////////////////////////////////////////////////////

// METReconstruction includes
#include "METReconstruction/METAssociator.h"
#include "xAODMissingET/MissingETComposition.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAssociationMap.h"
#include "AthContainers/ConstDataVector.h"

// Tracking EDM
#include "xAODTracking/Vertex.h"
#include "xAODTracking/TrackParticle.h"

// Track errors
#include "EventPrimitives/EventPrimitivesHelpers.h"

// Tool interface headers
#include "InDetTrackSelectionTool/IInDetTrackSelectionTool.h"
#include "RecoToolInterfaces/ITrackIsolationTool.h"
#include "RecoToolInterfaces/ICaloTopoClusterIsolationTool.h"

// For DeltaR
#include "FourMomUtils/xAODP4Helpers.h"

namespace met {

  using namespace xAOD;

  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

  // Constructors
  ////////////////
  METAssociator::METAssociator(const std::string& name) :
    AsgTool(name),
    m_trkseltool(this,""),
    m_trkIsolationTool(this,""),
    m_caloIsolationTool(this,"")
  {
    ATH_MSG_INFO("METAssoc constructor");
    declareProperty( "UseModifiedClus",    m_useModifiedClus = false             );
    declareProperty( "UseTracks",          m_useTracks   = true                  );
    declareProperty( "PFlow",              m_pflow       = false                 );
    declareProperty( "UseRapidity",        m_useRapidity = false                 );
    declareProperty( "TrackSelectorTool",  m_trkseltool                          );
    declareProperty( "TrackIsolationTool", m_trkIsolationTool                    );
    declareProperty( "CaloIsolationTool",  m_caloIsolationTool                   );
    declareProperty( "IgnoreJetConst",     m_skipconst = false                   );
    declareProperty( "ForwardColl",        m_forcoll   = ""                      );
    declareProperty( "ForwardDef",         m_foreta    = 2.5                     );
    declareProperty( "CentralTrackPtThr",  m_cenTrackPtThr = 30e+3               );
    declareProperty( "ForwardTrackPtThr",  m_forTrackPtThr = 30e+3               );
    declareProperty( "CleanCPFO",          m_cleanChargedPFO = true              );
    declareProperty( "UsePFOLinks",        m_usePFOLinks = false                 ); 
    declareProperty( "UseFELinks",         m_useFELinks = false                  ); 
    declareProperty( "NeutralPFOLinksKey", m_neutralPFOLinksKey = "neutralpfoLinks"); 
    declareProperty( "ChargedPFOLinksKey", m_chargedPFOLinksKey = "chargedpfoLinks"); 
    declareProperty( "NeutralFELinksKey",  m_neutralFELinksKey  = "neutralFELinks"); 
    declareProperty( "ChargedFELinksKey",  m_chargedFELinksKey  = "chargedFELinks"); 
  }

  // Destructor
  ///////////////
  METAssociator::~METAssociator()
  {}

  // Athena algtool's Hooks
  ////////////////////////////
  StatusCode METAssociator::initialize()
  {
    ATH_MSG_DEBUG ("Initializing " << name() << "...");

    ATH_CHECK( m_trkseltool.retrieve() );
    ATH_CHECK(m_trkIsolationTool.retrieve());
    ATH_CHECK(m_caloIsolationTool.retrieve());

    if(m_clcollKey.key() == "CaloCalTopoClusters") {
      if(m_useModifiedClus) {
        ATH_MSG_ERROR("Configured to use modified topocluster collection but \"CaloCalTopoClusters\" collection specified!");
        return StatusCode::FAILURE;
      } else {
        ATH_MSG_INFO("Configured to use standard topocluster collection.");
      }
    } else {
      if(m_useModifiedClus) {
        ATH_MSG_INFO("Configured to use modified topocluster collection \"" << m_clcollKey.key() << "\".");
      } else {
        ATH_MSG_ERROR("Configured to use topocluster collection \"" << m_clcollKey.key() << "\", but modified clusters flag not set!");
        return StatusCode::FAILURE;
      }
    }

    //initialise read handle keys
    ATH_CHECK( m_pvcollKey.initialize(m_useTracks));
    ATH_CHECK( m_trkcollKey.initialize(m_useTracks));
    
    ATH_CHECK(m_fecollKey.initialize(m_pflow && !m_fecollKey.key().empty()));
    ATH_CHECK( m_pfcollKey.initialize(m_pflow && m_fecollKey.key().empty()));
    if(m_pflow){
      if(!m_fecollKey.key().empty()){
        ATH_MSG_INFO("Configured to use FlowElement collection \"" << m_fecollKey.key() << "\".");
      }
      else{
        ATH_MSG_INFO("Configured to use PFlow collection \"" << m_pfcollKey.key() << "\".");
        if(m_pfcollKey.key() == "JetETMissParticleFlowObjects") {
          ATH_MSG_ERROR("Configured to use standard pflow collection \"" << m_pfcollKey.key() << "\".");
          ATH_MSG_ERROR("This is no longer supported -- please use the CHSParticleFlowObjects collection, which has the four-vector corrections built in.");
          return StatusCode::FAILURE;
        }
      }
    }

    ATH_CHECK( m_clcollKey.initialize(!m_skipconst || m_forcoll.empty()));

    std::string hybridname = "Etmiss";
    hybridname += m_clcollKey.key();
    hybridname += m_foreta;
    hybridname += m_forcoll;
    ATH_CHECK( m_hybridContKey.assign(hybridname));
    ATH_CHECK( m_hybridContKey.initialize(m_skipconst && !m_forcoll.empty()));

    return StatusCode::SUCCESS;
  }

  StatusCode METAssociator::execute(xAOD::MissingETContainer* metCont, xAOD::MissingETAssociationMap* metMap) const
  {
    ATH_MSG_DEBUG ("In execute: " << name() << "...");
    if(!metCont) {
      ATH_MSG_WARNING("Invalid pointer to MissingETContainer supplied! Abort.");
      return StatusCode::FAILURE;
    }

    if(!metMap) {
      ATH_MSG_WARNING("Invalid pointer to MissingETAssociationMap supplied! Abort.");
      return StatusCode::FAILURE;
    }
    if(m_pflow && !m_useTracks ){
      ATH_MSG_WARNING("Attempting to build PFlow MET without a track collection.");
      return StatusCode::FAILURE;
    }

    return this->executeTool(metCont, metMap);
  }

  StatusCode METAssociator::retrieveConstituents(met::METAssociator::ConstitHolder& constits) const
  {
    ATH_MSG_DEBUG ("In execute: " << name() << "...");
    if (!m_skipconst || m_forcoll.empty()) {

      SG::ReadHandle<IParticleContainer> topoclusterCont(m_clcollKey);
      if (!topoclusterCont.isValid()) {
        ATH_MSG_WARNING("Unable to retrieve topocluster container " << m_clcollKey.key() << " for overlap removal");
        return StatusCode::FAILURE;
      }
      constits.tcCont=topoclusterCont.cptr();
      ATH_MSG_DEBUG("Successfully retrieved topocluster collection");
    } else {
      std::string hybridname = "Etmiss";
      hybridname += m_clcollKey.key();
      hybridname += m_foreta;
      hybridname += m_forcoll;

      SG::ReadHandle<IParticleContainer> hybridCont(m_hybridContKey);
      if( hybridCont.isValid()) {
        constits.tcCont=hybridCont.cptr();
      } else {
        ATH_MSG_WARNING("Trying to do something currently unsupported- lets abort");
        return StatusCode::FAILURE;
        // Trying to do this using write handles (need to get some input here)
        /*std::unique_ptr<ConstDataVector<IParticleContainer>> hybridCont = std::make_unique<ConstDataVector<IParticleContainer>>();
        SG::WriteHandle<ConstDataVector<IParticleContainer>> hybridContHandle(hybridname);

        StatusCode sc = hybridContHandle.record(std::make_unique<ConstDataVector<IParticleContainer>>(*hybridCont));

        if (sc.isFailure()) {
          ATH_MSG_WARNING("Unable to record container");
         return StatusCode::SUCCESS;

        }*/

        /*SG::ReadHandle<IParticleContainer> centCont(m_clcoll);
        if (!centCont.isValid()) {
          ATH_MSG_WARNING("Unable to retrieve central container " << m_clcoll << " for overlap removal");
          return StatusCode::FAILURE;
        }

        SG::ReadHandle<IParticleContainer> forCont(m_forcoll);
        if (!forCont.isValid()) {
          ATH_MSG_WARNING("Unable to retrieve forward container " << m_forcoll << " for overlap removal");
          return StatusCode::FAILURE;
        }
        ConstDataVector<IParticleContainer> *hybridCont = new ConstDataVector<IParticleContainer>(SG::VIEW_ELEMENTS);

        const IParticleContainer* centCont=0;
        if( evtStore()->retrieve(centCont, m_clcoll).isFailure() ) {
          ATH_MSG_WARNING("Unable to retrieve central container " << m_clcoll << " for overlap removal");
          return StatusCode::FAILURE;
        }

        const IParticleContainer* forCont=0;
        if( evtStore()->retrieve(forCont, m_forcoll).isFailure() ) {
          ATH_MSG_WARNING("Unable to retrieve forward container " << m_forcoll << " for overlap removal");
          return StatusCode::FAILURE;
        }

        for(const auto clus : *centCont) if (fabs(clus->eta())<m_foreta) hybridCont->push_back(clus);
        for(const auto clus : *forCont) if (fabs(clus->eta())>=m_foreta) hybridCont->push_back(clus);
        ATH_CHECK( evtStore()->record(hybridCont,hybridname));
        constits.tcCont = hybridCont->asDataVector();
        */
      }
    }

    if( !m_useTracks){
      //if you want to skip tracks, set the track collection empty manually
      ATH_MSG_DEBUG("Skipping tracks");
    }else{
      SG::ReadHandle<VertexContainer> vxCont(m_pvcollKey);
      if (!vxCont.isValid()) {
        ATH_MSG_WARNING("Unable to retrieve primary vertex container " << m_pvcollKey.key());
        //this is actually really bad.  If it's empty that's okay
        return StatusCode::FAILURE;
      }

      ATH_MSG_DEBUG("Successfully retrieved primary vertex container");
      ATH_MSG_DEBUG("Container holds " << vxCont->size() << " vertices");

      for(const auto vx : *vxCont) {
        ATH_MSG_VERBOSE( "Testing vertex " << vx->index() );
        if(vx->vertexType()==VxType::PriVtx)
          {constits.pv = vx; break;}
      }
      if(!constits.pv) {
        ATH_MSG_DEBUG("Failed to find primary vertex! Reject all tracks.");
      } else {
        ATH_MSG_VERBOSE("Primary vertex has z = " << constits.pv->z());
      }

      constits.trkCont=0;
      ATH_MSG_DEBUG("Retrieving Track collection " << m_trkcollKey.key());
      SG::ReadHandle<TrackParticleContainer> trCont(m_trkcollKey);
      if (!trCont.isValid()) {
        ATH_MSG_WARNING("Unable to retrieve track particle container");
        return StatusCode::FAILURE;
      }
      constits.trkCont=trCont.cptr();

      if(m_pflow){
        if(!m_fecollKey.key().empty()){
          ATH_MSG_DEBUG("Retrieving FlowElement collection " << m_fecollKey.key());
          constits.feCont = 0;
          SG::ReadHandle<xAOD::FlowElementContainer> feCont(m_fecollKey);
          if (!feCont.isValid()) {
            ATH_MSG_ERROR("Unable to retrieve FlowElement container "<< m_fecollKey.key());
            return StatusCode::FAILURE;
          }
          constits.feCont=feCont.cptr();
        }
        else{
          ATH_MSG_DEBUG("Retrieving PFlow collection " << m_pfcollKey.key());
          constits.pfoCont = 0;
          SG::ReadHandle<PFOContainer> pfCont(m_pfcollKey);
          if (!pfCont.isValid()) {
            ATH_MSG_WARNING("Unable to PFlow object container");
            return StatusCode::FAILURE;
          }
          constits.pfoCont=pfCont.cptr();
        }
      }//pflow
    }//retrieve track/pfo containers

    return StatusCode::SUCCESS;
  }

  ///////////////////////////////////////////////////////////////////
  // Protected methods:
  ///////////////////////////////////////////////////////////////////

  StatusCode METAssociator::fillAssocMap(xAOD::MissingETAssociationMap* metMap,
                                         const xAOD::IParticleContainer* hardObjs) const
  {
    ConstitHolder constits;

    if (retrieveConstituents(constits).isFailure()) {
      ATH_MSG_DEBUG("Unable to retrieve constituent containers");
      return StatusCode::FAILURE;
    }

    std::vector<const IParticle*> constlist;
    constlist.reserve(20);
    std::vector<const IParticle*> hardObjs_tmp;
    for(const auto obj : *hardObjs) {
      hardObjs_tmp.push_back(obj);
    }
    std::sort(hardObjs_tmp.begin(),hardObjs_tmp.end(),greaterPt);

    for(const auto& obj : hardObjs_tmp) {
      if(obj->pt()<5e3 && obj->type()!=xAOD::Type::Muon) continue;
      constlist.clear();
      ATH_MSG_VERBOSE( "Object type, pt, eta, phi = " << obj->type() << ", " << obj->pt() << ", " << obj->eta() << "," << obj->phi() );
      if(m_pflow){
        if(!m_fecollKey.key().empty()){
          if(!m_useTracks){
            ATH_MSG_ERROR("Attempting to build FlowElement MET without a track collection.");
            return StatusCode::FAILURE;
          }
          std::map<const IParticle*, MissingETBase::Types::constvec_t> momentumOverride;
          ATH_CHECK( this->extractFE(obj, constlist, constits, momentumOverride) );
          MissingETComposition::insert(metMap, obj, constlist, momentumOverride);
        }
        else{
          // Old PFO EDM
          if(!m_useTracks){
            ATH_MSG_DEBUG("Attempting to build PFlow without a track collection.");
            return StatusCode::FAILURE;
          }else{
            std::map<const IParticle*,MissingETBase::Types::constvec_t> momentumOverride;
            ATH_CHECK( this->extractPFO(obj,constlist,constits,momentumOverride) );
            MissingETComposition::insert(metMap,obj,constlist,momentumOverride);
          }
        }
      } else {
        std::vector<const IParticle*> tclist;
        tclist.reserve(20);
        ATH_CHECK( this->extractTopoClusters(obj,tclist,constits) );
        if(m_useModifiedClus) {
          for(const auto& cl : tclist) {
            // use index-parallelism to identify shallow copied constituents
            constlist.push_back((*constits.tcCont)[cl->index()]);
          }
        } else {
          constlist = tclist;
        }
        if(m_useTracks) ATH_CHECK( this->extractTracks(obj,constlist,constits) );
        MissingETComposition::insert(metMap,obj,constlist);
      }
    }
    return StatusCode::SUCCESS;
  }

  // Accept Track
  ////////////////
  bool METAssociator::acceptTrack(const xAOD::TrackParticle* trk, const xAOD::Vertex* vx) const
  {

    if (!vx) return false;//in events with no pv, we will just reject all tracks, and therefore build only the calo MET
    return static_cast<bool>  (m_trkseltool->accept( *trk, vx ));
 }


  bool METAssociator::isGoodEoverP(const xAOD::TrackParticle* trk) const
  {

    if( (fabs(trk->eta())<1.5 && trk->pt()>m_cenTrackPtThr) ||
        (fabs(trk->eta())>=1.5 && trk->pt()>m_forTrackPtThr) ) {

      // Get relative error on qoverp
      float Rerr = Amg::error(trk->definingParametersCovMatrix(),4)/fabs(trk->qOverP());
      ATH_MSG_VERBOSE( "Track momentum error (%): " << Rerr*100 );

      // first compute track and calo isolation variables
      float ptcone20 = 0., isolfrac = 0., etcone10 = 0., EoverP = 0.;
      // ptcone
      TrackIsolation trkIsoResult;
      std::vector<Iso::IsolationType> trkIsoCones; 
      trkIsoCones.push_back(xAOD::Iso::IsolationType::ptcone20);
      xAOD::TrackCorrection trkIsoCorr;
      trkIsoCorr.trackbitset.set(xAOD::Iso::IsolationTrackCorrection::coreTrackPtr); 
      m_trkIsolationTool->trackIsolation(trkIsoResult,
                                         *trk,
                                         trkIsoCones,
                                         trkIsoCorr);
      ptcone20 = trkIsoResult.ptcones.size() > 0 ? trkIsoResult.ptcones[0] : 0;
      isolfrac = ptcone20/trk->pt();
      // etcone
      CaloIsolation caloIsoResult;
      std::vector<Iso::IsolationType> caloIsoCones; 
      // We can't actually configure the tool to give etcone10, so instead we have to compute etcone20,
      // applying the core cone correction.
      // Then, we retrieve the correction value, which is etcone10, rather than the isolation value
      caloIsoCones.push_back(xAOD::Iso::IsolationType::etcone20); 
      xAOD::CaloCorrection caloIsoCorr_coreCone;
      caloIsoCorr_coreCone.calobitset.set(xAOD::Iso::IsolationCaloCorrection::coreCone); // this is etcone10
      m_caloIsolationTool->caloTopoClusterIsolation(caloIsoResult,
                                                    *trk,
                                                    caloIsoCones,
                                                    caloIsoCorr_coreCone);
      if(caloIsoResult.etcones.size() > 0) {
        // retrieve the correction value for the core cone
        etcone10 = caloIsoResult.coreCorrections[xAOD::Iso::IsolationCaloCorrection::coreCone][xAOD::Iso::IsolationCorrectionParameter::coreEnergy];
      } else {
        ATH_MSG_WARNING("isGoodEoverP: Failed to retrieve the isolation core correction (etcone10)! Setting etcone10=0");
        etcone10 = 0.;
      }
      EoverP   =  etcone10/trk->pt(); 
      /////////////////////////////////////////////////////////////////////////
      ATH_MSG_VERBOSE( "Track isolation fraction: " << isolfrac );
      ATH_MSG_VERBOSE( "Track E/P = " << EoverP );

      if(isolfrac<0.1) {
            // isolated track cuts
            if(Rerr>0.4) return false;
            else if (EoverP<0.65 && ((EoverP>0.1 && Rerr>0.05) || Rerr>0.1)) return false;
      } else {
            // non-isolated track cuts
            float trkptsum = ptcone20+trk->pt();
            if(etcone10/trkptsum<0.6 && trk->pt()/trkptsum>0.6) return false;
      }
    }
    return true;
  }

}
