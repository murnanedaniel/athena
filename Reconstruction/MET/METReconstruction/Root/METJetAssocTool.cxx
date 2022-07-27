///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// METJetAssocTool.cxx
// Implementation file for class METJetAssocTool
//
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//
// Author: P Loch, S Resconi, TJ Khoo, AS Mete
///////////////////////////////////////////////////////////////////

// METReconstruction includes
#include "METReconstruction/METJetAssocTool.h"

// MET EDM
#include "xAODMissingET/MissingETComposition.h"

// Jet EDM
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAttributes.h"
#include "xAODJet/JetTypes.h"

// Tracking EDM
#include "xAODTracking/Vertex.h"

// Cluster EDM
#include "xAODCaloEvent/CaloCluster.h"

// DeltaR calculation
#include "FourMomUtils/xAODP4Helpers.h"

namespace met {

  using namespace xAOD;

  // Constructors
  ////////////////
  METJetAssocTool::METJetAssocTool(const std::string& name) :
    AsgTool(name),
    METAssociator(name)
  {
  }

  // Athena algtool's Hooks
  ////////////////////////////
  StatusCode METJetAssocTool::initialize()
  {
    ATH_CHECK( METAssociator::initialize() );
    ATH_MSG_VERBOSE ("Initializing " << name() << "...");
    //Initialise ReadHandles
    ATH_CHECK( m_jetContKey.initialize());


    return StatusCode::SUCCESS;
  }

  StatusCode METJetAssocTool::executeTool(xAOD::MissingETContainer* /*metCont*/, xAOD::MissingETAssociationMap* metMap) const
  {
    ATH_MSG_VERBOSE ("In execute: " << name() << "...");

    // Retrieve the jet container
    SG::ReadHandle<xAOD::JetContainer> jetCont(m_jetContKey);
    if (!jetCont.isValid()) {
      ATH_MSG_WARNING("Unable to retrieve input jet container " << m_jetContKey.key());
      return StatusCode::FAILURE;
    }

    ATH_MSG_DEBUG("Successfully retrieved jet collection");

    ConstitHolder constits;
    if (retrieveConstituents(constits).isFailure()) {
      ATH_MSG_WARNING("Unable to retrieve constituent containers");
      return StatusCode::FAILURE;
    }

    std::set<const xAOD::IParticle*> newConst;
    if (m_skipconst) {
      for (const auto clust : *constits.tcCont) newConst.insert(clust);
    }

    // Create jet associations
    for(const auto jet : *jetCont) {
      std::vector<const IParticle*> selectedTracks;
      bool mismatchedPFlow = m_pflow && (jet->rawConstituent(0)->type()!=xAOD::Type::ParticleFlow && jet->rawConstituent(0)->type()!=xAOD::Type::FlowElement);
      bool mismatchedState = !m_skipconst && !m_pflow && jet->rawConstituent(0)->type()==xAOD::Type::CaloCluster && ((static_cast<const xAOD::CaloCluster*>(jet->rawConstituent(0))->signalState()==xAOD::CaloCluster::CALIBRATED && jet->getConstituentsSignalState()==xAOD::UncalibratedJetConstituent) || (static_cast<const xAOD::CaloCluster*>(jet->rawConstituent(0))->signalState()==xAOD::CaloCluster::UNCALIBRATED && jet->getConstituentsSignalState()==xAOD::CalibratedJetConstituent));
      bool newConstVec = m_skipconst || mismatchedPFlow || mismatchedState;
      if (m_pflow && !mismatchedPFlow) {
        if(jet->rawConstituent(0)->type() == xAOD::Type::FlowElement){
          for (size_t consti = 0; consti < jet->numConstituents(); consti++) {
            const xAOD::FlowElement *pfo = static_cast<const xAOD::FlowElement*>(jet->rawConstituent(consti));
            ATH_MSG_VERBOSE("Jet constituent PFO, pt :" << pfo->pt() << ", charge: " << pfo->charge());
            if ((pfo->isCharged()) && (!m_cleanChargedPFO || isGoodEoverP(static_cast<const xAOD::TrackParticle*>(pfo->chargedObject(0))))) {
              ATH_MSG_VERBOSE("  Accepted charged PFO, pt " << pfo->pt());
              selectedTracks.push_back(pfo);
            }
          }
        }
        else{
          // constituents are xAOD::PFO
          for (size_t consti = 0; consti < jet->numConstituents(); consti++) {
            const xAOD::PFO *pfo = static_cast<const xAOD::PFO*>(jet->rawConstituent(consti));
            ATH_MSG_VERBOSE("Jet constituent PFO, pt :" << pfo->pt() << ", charge: " << pfo->charge());
            if (pfo->isCharged() && (!m_cleanChargedPFO || isGoodEoverP(pfo->track(0)))) {
              ATH_MSG_VERBOSE("  Accepted charged PFO, pt " << pfo->pt());
              selectedTracks.push_back(pfo);
            }
          }
        }
      } else {
        std::vector<const IParticle*> jettracks;
        jet->getAssociatedObjects<IParticle>(JetAttribute::GhostTrack,jettracks);

        selectedTracks.reserve(jettracks.size());
        for(const auto& trk : jettracks) {
          const TrackParticle* pTrk = static_cast<const TrackParticle*>(trk);
          if(acceptTrack(pTrk,constits.pv) && isGoodEoverP(pTrk)) {
            selectedTracks.push_back(trk);
            ATH_MSG_VERBOSE("Accept track " << trk << " px, py = " << trk->p4().Px() << ", " << trk->p4().Py());
          }
        }
      }
      std::vector<const IParticle*> consts;
      std::map<const IParticle*,MissingETBase::Types::constvec_t> momenta;

      MissingETComposition::add(metMap,jet,selectedTracks);
      if (mismatchedPFlow) getPFOs(jet,consts,constits,momenta);
      else if (mismatchedState) getClus(jet,consts);
      else if (newConstVec) getOther(jet,consts,&newConst);
      if (newConstVec) MissingETComposition::setJetConstSum(metMap,jet,consts,momenta);
      ATH_MSG_VERBOSE("Added association " << metMap->findIndex(jet) << " pointing to jet " << jet);
      ATH_MSG_VERBOSE("Jet pt, eta, phi = " << jet->pt() << ", " << jet->eta() << "," << jet->phi() );

    }
    MissingETComposition::addMiscAssociation(metMap);
    ATH_MSG_DEBUG("Added miscellaneous association");

    return StatusCode::SUCCESS;
  }

  void METJetAssocTool::getPFOs(const xAOD::Jet *jet,
                                std::vector<const xAOD::IParticle*> &consts,
                                const met::METAssociator::ConstitHolder& constits,
                                std::map<const xAOD::IParticle*,MissingETBase::Types::constvec_t> &momenta) const {

    std::vector<const IParticle*> jettracks;
    jet->getAssociatedObjects<IParticle>(JetAttribute::GhostTrack,jettracks);

    if(constits.feCont != nullptr){
      // We have PFOs in FlowElement format
      for(const xAOD::FlowElement* pfo : *constits.feCont) {
        if (pfo->isCharged()) {
          const TrackParticle* pfotrk = static_cast<const xAOD::TrackParticle*>(pfo->chargedObject(0));
          for(const xAOD::IParticle* track : jettracks) {
            const TrackParticle* jettrk = static_cast<const xAOD::TrackParticle*>(track);
            if (jettrk==pfotrk) {
              consts.push_back(pfo);
              break;
            }
          }
        }
        else {
          bool marked = false;
          for (size_t consti = 0; consti < jet->numConstituents(); consti++){
            if (pfo->p4().DeltaR(jet->rawConstituent(consti)->p4())<0.05) marked = true;
          }
          if (marked) {
            consts.push_back(pfo);
            TLorentzVector momentum = pfo->p4();
            momenta[pfo] = MissingETBase::Types::constvec_t(momentum.Px(),momentum.Py(),momentum.Pz(),
                                                            momentum.E(),momentum.Pt());
          }
        }
      }
    }
    else{
      // No FlowElements, assume xAOD::PFO format
      for(const auto pfo : *constits.pfoCont) {
        if (pfo->isCharged()) {
          const TrackParticle* pfotrk = pfo->track(0);
          for(const auto& trk : jettracks) {
            if (trk==pfotrk) {
              consts.push_back(pfo);
              break;
            }
          }
        } else {
          bool marked = false;
          for (size_t consti = 0; consti < jet->numConstituents(); consti++) if (pfo->p4().DeltaR(jet->rawConstituent(consti)->p4())<0.05) marked = true;
          if (marked) {
            consts.push_back(pfo);
            TLorentzVector momentum = pfo->p4();
            momenta[pfo] = MissingETBase::Types::constvec_t(momentum.Px(),momentum.Py(),momentum.Pz(),
                                                            momentum.E(),momentum.Pt());
          }
        }
      }
    }
  }

  void METJetAssocTool::getClus(const xAOD::Jet *jet,
                                std::vector<const xAOD::IParticle*> &consts) const {
    std::vector<ElementLink<IParticleContainer> > jetconst = jet->constituentLinks();
    for(const auto& clus : jetconst) consts.push_back(*clus);
  }

  void METJetAssocTool::getOther(const xAOD::Jet *jet,
                                 std::vector<const xAOD::IParticle*> &consts,
                                 std::set<const xAOD::IParticle*> *newConst) const {
    std::vector<ElementLink<IParticleContainer> > jetconst = jet->constituentLinks();
    for(const auto& clus : *newConst) if (clus->container()!=jet->rawConstituent(0)->container() && clus->e()>0 && xAOD::P4Helpers::isInDeltaR(*jet,*clus,m_matchRadius,m_useRapidity)) consts.push_back(clus);
    //for(const auto& clus : *newConst) if (clus->type()!=jet->rawConstituent(0)->type() && clus->e()>0 && xAOD::P4Helpers::isInDeltaR(*jet,*clus,m_matchRadius,m_useRapidity)) consts.push_back(clus);
    for(const auto& clusL : jetconst) {
      const xAOD::IParticle *clus = *clusL;
      if (newConst->count(clus)) consts.push_back(clus);
    }
  }
}
