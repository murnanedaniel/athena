////////////////////-*- C++ -*-////////////////////////////////////

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

// PFlowAugmentationTool.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef DERIVATIONFRAMEWORK_PFLOWAUGMENTATIONTOOL_H
#define DERIVATIONFRAMEWORK_PFLOWAUGMENTATIONTOOL_H

#include <string>
#include <vector>

#include "AthenaBaseComps/AthAlgTool.h"
#include "DerivationFrameworkInterfaces/IAugmentationTool.h"
#include "GaudiKernel/ToolHandle.h"

#include "StoreGate/WriteDecorHandleKey.h"

#include "xAODTracking/Vertex.h"
#include "xAODTracking/VertexContainer.h"
#include <xAODPFlow/FlowElementContainer.h>
#include <PFlowUtils/IWeightPFOTool.h>

namespace DerivationFramework {

  class PFlowAugmentationTool : public AthAlgTool, public IAugmentationTool {
  public: 
    PFlowAugmentationTool(const std::string& t, const std::string& n, const IInterface* p);

    StatusCode initialize();
    StatusCode finalize();
    virtual StatusCode addBranches() const;

  private:

    float m_z0sinthcut;

    ToolHandle<CP::IWeightPFOTool> m_weightPFOTool;    /// Retrieval tool

    SG::ReadHandleKey<xAOD::VertexContainer> m_vertexContainer_key{this, "VertexContainer", "PrimaryVertices", "Input vertex container"};
    SG::ReadHandleKey<xAOD::FlowElementContainer> m_pfoContainer_key{this, "GlobalChargedParticleFlowObjects", "GlobalChargedParticleFlowObjects", "Input charged PFO"};

    SG::WriteDecorHandleKey<xAOD::FlowElementContainer> m_PVmatchedKey {this, "PVmatchedKey", "GlobalChargedParticleFlowObjects.DFCommonPFlow_PVMatched", "Boolean indicating if PFO was matched to PV "};
    SG::WriteDecorHandleKey<xAOD::FlowElementContainer> m_corrP4_ptKey {this, "m_corrP4_ptKey", "GlobalChargedParticleFlowObjects.DFCommonPFlow_CaloCorrectedPt", "Decoration for weighted charged PFO pt"};
    SG::WriteDecorHandleKey<xAOD::FlowElementContainer> m_z0Key {this, "m_z0Key", "GlobalChargedParticleFlowObjects.DFCommonPFlow_z0", "Decoration for track z0"};
    SG::WriteDecorHandleKey<xAOD::FlowElementContainer> m_vzKey{this, "m_vzKey","GlobalChargedParticleFlowObjects.DFCommonPFlow_vz", "Decoration for track vz"};
    SG::WriteDecorHandleKey<xAOD::FlowElementContainer> m_d0Key{this, "m_d0Key","GlobalChargedParticleFlowObjects.DFCommonPFlow_d0", "Decoration for track d0"};
    SG::WriteDecorHandleKey<xAOD::FlowElementContainer> m_thetaKey{this, "m_thetaKey","GlobalChargedParticleFlowObjects.DFCommonPFlow_theta", "Decoration for track theta"};
    SG::WriteDecorHandleKey<xAOD::FlowElementContainer> m_envWeightKey{this, "m_envWeightKey","GlobalChargedParticleFlowObjects.DFCommonPFlow_envWeight", "Decoration for weight for dense environments"};
    
  }; 
}

#endif // DERIVATIONFRAMEWORK_PFLOWAUGMENTATIONTOOL_H
