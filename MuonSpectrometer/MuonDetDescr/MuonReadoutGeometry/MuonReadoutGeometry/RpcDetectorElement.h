/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/***************************************************************************
 Collect RPC readout elements belonging to the same DoubletR - granularity is same as for EDM (hash ids)
 -------------------------------------------------------
***************************************************************************/

#ifndef MUONREADOUTGEOMETRY_RPCDETECTORELEMENT_H
#define MUONREADOUTGEOMETRY_RPCDETECTORELEMENT_H

#include "MuonReadoutGeometry/MuonDetectorElement.h"

#include "Identifier/Identifier.h"
#include "Identifier/IdentifierHash.h"

#include <vector>

namespace Trk{
  class Surface;
  class SurfaceBounds;
}

class GeoVFullPhysVol;
class RpcIdHelper;

namespace MuonGM {
    
  class MuonDetectorManager;
  class MuonStation;
  class RpcReadoutElement;

  typedef std::vector<const RpcReadoutElement *> REVector;
  typedef std::vector<const RpcReadoutElement *>::const_iterator REIterator;

    
  class RpcDetectorElement: public MuonDetectorElement
  {

  public:

    RpcDetectorElement(GeoVFullPhysVol* pv, MuonDetectorManager* mgr, Identifier id, IdentifierHash idHash);
   
    virtual int getStationEta() const {return 0;}; //!< returns stationEta 
    virtual int getStationPhi() const {return 0;}; //!< returns stationPhi

    // access to the readout-elements in this DetectorElement
    const RpcReadoutElement* getRpcReadoutElement(Identifier id) const;
    // this is a channelId

    const RpcReadoutElement* getRpcReadoutElement(int dbz, int dbp) const;
    // access to the MuonStation this DetectorElement belongs to
    MuonStation* parentMuonStation() const;

    unsigned int NdoubletZ() const;
    unsigned int NsegmentedDoubletZ() const;
    unsigned int NPhimodules(int dbz) const;
    void addRpcReadoutElement(const RpcReadoutElement* rpc, int index);

    unsigned int nMDTinStation() const {return 0;} 
    unsigned int nCSCinStation() const {return 0;}
    unsigned int nTGCinStation() const {return 0;}
    unsigned int nRPCinStation() const {return nReadoutElements();}

    const Amg::Transform3D& transform() const;

    const Trk::Surface& surface() const;
  
    const Trk::SurfaceBounds& bounds() const;

    const Amg::Vector3D& center() const;
  
    const Amg::Vector3D& normal() const;
  
    const Amg::Vector3D& normal(const Identifier& id) const;
  
    const Trk::Surface& surface(const Identifier& id) const;
  
    const Trk::SurfaceBounds& bounds(const Identifier& id) const;
  
    const Amg::Transform3D& transform(const Identifier& id) const;
  
    const Amg::Vector3D& center(const Identifier& id) const;

    std::vector<const Trk::Surface*> surfaces() const;

    enum RpcGMRanges
      {NDoubletZ = 4}; 
    // using some trick to save space: dbz=4 if rib's chambers and doubletphi=2;

  private:

    const RpcIdHelper* m_helper;
    int m_ndbz;
    const RpcReadoutElement* m_rpcVector[NDoubletZ];
  };

} // namespace MuonGM

#endif // MUONREADOUTGEOMETRY_RPCDETECTORELEMENT_H
