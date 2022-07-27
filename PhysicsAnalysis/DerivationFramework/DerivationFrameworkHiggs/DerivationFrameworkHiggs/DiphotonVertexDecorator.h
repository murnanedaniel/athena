/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

//////////////////////////////////////////////////////////////////
// DiphotonVertexDecorator.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef DERIVATIONFRAMEWORK_DiphotonVertexDecorator_H
#define DERIVATIONFRAMEWORK_DiphotonVertexDecorator_H
 
#include <string>
#include <vector>
#include <algorithm>

// Gaudi & Athena basics
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"

// DerivationFramework includes
#include "DerivationFrameworkInterfaces/IAugmentationTool.h"
// xAOD header files
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTracking/VertexContainer.h"

#include "xAODPFlow/FlowElementContainer.h"

#include "PhotonVertexSelection/IPhotonVertexSelectionTool.h"

namespace DerivationFramework {

  /** @class DiphotonVertexDecorator
      @author Bruno Lenzi
      @author Leo Cerda
      @author magdac@cern.ch
     */
  class DiphotonVertexDecorator : public AthAlgTool, public IAugmentationTool {

    public: 
      /** Constructor with parameters */
      DiphotonVertexDecorator( const std::string& t, const std::string& n, const IInterface* p );

      /** Destructor */
      ~DiphotonVertexDecorator();

      // Athena algtool's Hooks
      StatusCode  initialize();
      StatusCode  finalize();

      virtual StatusCode addBranches() const;

    private:

      ///////////////
      ///// TOOLS 
      ToolHandle<CP::IPhotonVertexSelectionTool> m_photonVertexSelectionTool{this,  "PhotonVertexSelectionTool", "", ""};
      
      SG::ReadHandleKey<xAOD::VertexContainer> m_primaryVertexKey{this, "PrimaryVertexName", "PrimaryVertices", "" };
      SG::ReadHandleKey<xAOD::PhotonContainer> m_photonKey { this, "PhotonKey", "Photons", "" };
      SG::WriteHandleKey<xAOD::VertexContainer> m_diphotonVertexKey{this, "DiphotonVertexName", "HggPrimaryVertices", "" };
      SG::ReadHandleKey<xAOD::FlowElementContainer> m_FEContainerHandleKey{this,"PFOContainerName","JetETMissChargedParticleFlowObjects","ReadHandleKey for the PFO container"};
      ///////////////
      ///// SETTINGS

      double m_minPhotonPt;
      bool   m_removeCrack;
      double m_maxEta;
      bool   m_ignoreConv;
      std::string m_pfoToolName;
      double m_tcMatch_dR;
      double m_tcMatch_maxRat;

      bool  PhotonPreselect(const xAOD::Photon *ph) const;
      StatusCode matchFlowElement(const xAOD::Photon* eg,const xAOD::FlowElementContainer *pfoCont) const;
      static inline bool greaterPtFlowElement(const xAOD::FlowElement* part1, const xAOD::FlowElement* part2) {
        if (part1->charge()==0 && part2->charge()!=0) return false;
        if (part1->charge()!=0 && part2->charge()==0) return true;
        return part1->pt()>part2->pt();
    }

  }; 

}

#endif // DERIVATIONFRAMEWORK_DiphotonVertexDecorator_H
