/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DERIVATIONFRAMEWORK_PVTOOLS_H
#define DERIVATIONFRAMEWORK_PVTOOLS_H

#include "GaudiKernel/ToolHandle.h"
#include "xAODBPhys/BPhysHelper.h"
#include <vector>
#include "xAODEventInfo/EventInfo.h"
  
// Author: Adam Barton <abarton@SPAMMENOTTtttcern.ch>
namespace Trk {
  class V0Tools;
}
namespace Analysis{
  class PrimaryVertexRefitter;
}


namespace DerivationFramework {
  
  class BPhysPVTools {

  private:
       const Trk::V0Tools *m_v0Tools;
       const xAOD::EventInfo *m_EvtData;

       /// minimum number of tracks required in PVs considered
       size_t m_PV_minNTracks;

       bool m_3dCalc;
       
  public:
  
       BPhysPVTools(const Trk::V0Tools *v0Tools);
       BPhysPVTools(const Trk::V0Tools *v0Tools, const xAOD::EventInfo*);
       void SetSave3d(bool v) { m_3dCalc =v; }
       StatusCode FillCandExistingVertices(xAOD::VertexContainer* vtxContainer, const xAOD::VertexContainer* pvContainer, int DoVertexType);
       
       static void FillBPhysHelperNULL(xAOD::BPhysHelper &vtx, const xAOD::VertexContainer* PvContainer,
           xAOD::BPhysHelper::pv_type pvtype, bool do3d = false);
       
       StatusCode FillCandwithRefittedVertices(xAOD::VertexContainer* vtxContainer, const xAOD::VertexContainer* pvContainer,xAOD::VertexContainer* refPvContainer, const Analysis::PrimaryVertexRefitter* , size_t in_PV_max, int DoVertexType);
       
       void DecorateWithNULL(xAOD::VertexContainer* vtxContainer,const xAOD::VertexContainer* pvContainer, int DoVertexType) const;
       
       void DecorateWithDummyVertex(xAOD::VertexContainer* vtxContainer, const xAOD::VertexContainer* pvContainer, const xAOD::Vertex* Dummy, int DoVertexType, bool SetOrignal) const;
       
       ///Fills the BPhysHelper object with the standard parameters
       void FillBPhysHelper(xAOD::BPhysHelper &vtx, const xAOD::Vertex* refPV,const xAOD::VertexContainer* refPvContainer,
                    xAOD::BPhysHelper::pv_type pvtype, int) const;
    
       ///Returns the index integer of the vertex with the lowest Z in relation to the given vertex
       size_t FindLowZIndex(const xAOD::BPhysHelper &Obj,
			    const std::vector<const xAOD::Vertex*> &PVlist,
			    const size_t PV_minNTracks=0) const;
       ///Returns the index integer of the vertex with the lowest A0 in relation to the given vertex
       size_t FindLowA0Index(const xAOD::BPhysHelper &Obj,
			     const std::vector<const xAOD::Vertex*> &PVlist,
			     const size_t PV_minNTracks=0) const;
       
       static size_t FindHighPtIndex(const std::vector<const xAOD::Vertex*> &PVlist);
       
       /// Static method call with
       /// DerivationFramework::BPhysDerHelpers::GetGoodPV
       /// Returns a std::vector containing only PVs of type 1 and 3 - HighPt
       /// and Pileup, which have at least PV_minNTracks tracks.
       static std::vector<const xAOD::Vertex*> GetGoodPV(const xAOD::VertexContainer* pvContainer);
       
       /// Set the minimum number of tracks required for primary vertices to be
       /// considered for primary vertex association to a secondary vertex.
       /// Note that this requirement will not be applied for finding
       /// the vertex with the highest pT sum (FindHighPtIndex()) since
       /// it would possibly exclude this vertex which has been marked
       /// earlier in the tool chain.
       void SetMinNTracksInPV(size_t PV_minNTracks);

       /// Get the current beamspot position either from cache or from
       /// BeamCondSvc.
       /// Before processing a new event, make sure to call
       /// GetBeamSpot();
       [[nodiscard]] Amg::Vector3D GetBeamSpot() const;

       /// Find the index for the PV with the lowest distance in z of
       /// the SV's DOCA point w.r.t. the beamline and the PV.
       size_t FindLowZ0BAIndex(xAOD::BPhysHelper &obj,
			       const std::vector<const xAOD::Vertex*> &PVlist,
			       const size_t PV_minNTracks=0) const;
       /// Calculate the distance along z axis between the PV and
       ///  SV's DOCA point w.r.t. the beamline.
       double DistInZtoDOCA(xAOD::BPhysHelper &obj,
			    const xAOD::Vertex* vertex) const;
       /// Point of DOCA w.r.t. the beamline backward extrapolated
       /// along the B candidate's momentum direction. 
       Amg::Vector3D DocaExtrapToBeamSpot(xAOD::BPhysHelper &obj) const;

       static void PrepareVertexLinks(xAOD::Vertex* theResult,
               const xAOD::TrackParticleContainer* importedTrackCollection);
  }; // class BPhysPVTools

} // namespace DerivationFramework


#endif // DERIVATIONFRAMEWORK_PVTOOLS_H
