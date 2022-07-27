///////////////////////// -*- C++ -*- /////////////////////////////

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// METAssociator.h
// Header file for class METAssociator
//
// This is the base class for tools that construct MET terms
// from other object collections.
//
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
//
// Author: P Loch, S Resconi, TJ Khoo, AS Mete
///////////////////////////////////////////////////////////////////
#ifndef METRECONSTRUCTION_METASSOCIATOR_H
#define METRECONSTRUCTION_METASSOCIATOR_H

// STL includes
#include <string>

// FrameWork includes
#include "AsgTools/AsgTool.h"
#include "AsgTools/ToolHandle.h"

// METRecoInterface includes
#include "METRecoInterface/IMETAssocToolBase.h"

#include "xAODJet/JetContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/Vertex.h"
#include "xAODPFlow/PFOContainer.h"
#include "xAODPFlow/PFO.h"
#include "xAODPFlow/FlowElementContainer.h"

namespace InDet {
  class IInDetTrackSelectionTool;
}
 
namespace xAOD {
  class ITrackIsolationTool;
  class ICaloTopoClusterIsolationTool;
}

namespace met {
  class METAssociator
    : virtual public asg::AsgTool,
      virtual public IMETAssocToolBase
  {
 
    ///////////////////////////////////////////////////////////////////
    // Public methods:
    ///////////////////////////////////////////////////////////////////
    public:
 
    struct ConstitHolder {
      const xAOD::TrackParticleContainer* trkCont = 0;
      // Use IParticleContainer for flexibility e.g. if combining clusters & towers
      const xAOD::IParticleContainer* tcCont = 0;
      const xAOD::PFOContainer* pfoCont = 0;
      const xAOD::FlowElementContainer* feCont = 0;
      const xAOD::Vertex* pv = 0;
    };

    // Constructor w/ name
    METAssociator(const std::string& name);
    // Default Destructor
    virtual ~METAssociator();

    // AsgTool Handles
    virtual StatusCode initialize() override;
    virtual StatusCode execute (xAOD::MissingETContainer* metCont, xAOD::MissingETAssociationMap* metMap) const override;

    ///////////////////////////////////////////////////////////////////
    // Protected methods:
    ///////////////////////////////////////////////////////////////////
    protected:

    ToolHandle<InDet::IInDetTrackSelectionTool> m_trkseltool;
    ToolHandle<xAOD::ITrackIsolationTool> m_trkIsolationTool;
    ToolHandle<xAOD::ICaloTopoClusterIsolationTool> m_caloIsolationTool;

    std::string m_neutralFELinksKey; 
    std::string m_chargedFELinksKey; 
    std::string m_neutralPFOLinksKey; 
    std::string m_chargedPFOLinksKey; 
    bool m_usePFOLinks; 
    bool m_useFELinks; 

    SG::ReadHandleKey<xAOD::VertexContainer>  m_pvcollKey{this,"PrimVxColl","PrimaryVertices","Primary Vertex Collection"};
    SG::ReadHandleKey<xAOD::IParticleContainer>  m_clcollKey{this,"ClusColl","CaloCalTopoClusters","Topo cluster Collection"};
    SG::ReadHandleKey<xAOD::TrackParticleContainer>  m_trkcollKey{this,"TrkColl","InDetTrackParticles","Track particle Collection"};
    SG::ReadHandleKey<xAOD::PFOContainer>  m_pfcollKey{this,"PFlowColl","","PFO Collection"};
    SG::ReadHandleKey<xAOD::FlowElementContainer>  m_fecollKey{this,"FlowElementCollection","","FlowElement Collection (overrides PFO if not empty)"};
    SG::ReadHandleKey<xAOD::IParticleContainer>  m_hybridContKey{this,"HybridKey","","Hybrid Collection"};

    bool m_pflow;
    bool m_useTracks;
    bool m_useRapidity;
    bool m_useIsolationTools;
    bool m_useModifiedClus;
    bool m_weight_charged_pfo;
    bool m_cleanChargedPFO;

    bool m_skipconst;
    std::string m_forcoll;
    double m_foreta;

    double m_cenTrackPtThr;
    double m_forTrackPtThr;



    // reconstruction process to be defined in the individual tools
    // pure virtual -- we have no default
    virtual StatusCode executeTool(xAOD::MissingETContainer* metCont, xAOD::MissingETAssociationMap* metMap) const = 0;
    StatusCode retrieveConstituents(met::METAssociator::ConstitHolder& constits) const;

    bool acceptTrack (const xAOD::TrackParticle* trk, const xAOD::Vertex* pv) const;
    bool isGoodEoverP(const xAOD::TrackParticle* trk) const;

    virtual StatusCode fillAssocMap(xAOD::MissingETAssociationMap* metMap,
                                    const xAOD::IParticleContainer* hardObjs) const;
    virtual StatusCode extractPFO(const xAOD::IParticle* obj,
                                  std::vector<const xAOD::IParticle*>& pfolist,
                                  const met::METAssociator::ConstitHolder& constits,
                                  std::map<const xAOD::IParticle*,MissingETBase::Types::constvec_t> &momenta) const = 0;
    virtual StatusCode extractFE(const xAOD::IParticle* obj,
                                 std::vector<const xAOD::IParticle*>& felist,
                                 const met::METAssociator::ConstitHolder& constits,
                                 std::map<const xAOD::IParticle*,MissingETBase::Types::constvec_t> &momenta) const = 0;
    virtual StatusCode extractTracks(const xAOD::IParticle* obj,
                                     std::vector<const xAOD::IParticle*>& constlist,
                                     const met::METAssociator::ConstitHolder& constits) const = 0;
    virtual StatusCode extractTopoClusters(const xAOD::IParticle* obj,
                                           std::vector<const xAOD::IParticle*>& tclist,
                                           const met::METAssociator::ConstitHolder& constits) const = 0;
    static inline bool greaterPt(const xAOD::IParticle* part1, const xAOD::IParticle* part2) {
      return part1->pt()>part2->pt();
    }
    static inline bool greaterPtPFO(const xAOD::PFO* part1, const xAOD::PFO* part2) {
      if (part1->charge()==0 && part2->charge()!=0) return false;
      if (part1->charge()!=0 && part2->charge()==0) return true;
      if (part1->charge()==0 && part2->charge()==0) return part1->ptEM()>part2->ptEM();
      return part1->pt()>part2->pt();
    }
    static inline bool greaterPtFE(const xAOD::FlowElement* part1, const xAOD::FlowElement* part2) {
      if (!(part1->isCharged()) && part2->isCharged()) return false;
      if (part1->isCharged() && !(part2->isCharged())) return true;
      return part1->pt() > part2->pt();
    }
    ///////////////////////////////////////////////////////////////////
    // Private methods:
    ///////////////////////////////////////////////////////////////////
    private:

    // Default Constructor
    METAssociator();
  };
}

#endif // METRECONSTRUCTION_METASSOCBUILDERTOOL_H
