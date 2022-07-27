/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// HGTD_LayerBuilderCond.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef HGTD_TRACKINGGEOMETRY_HGTDLAYERBUILDERCOND_H
#define HGTD_TRACKINGGEOMETRY_HGTDLAYERBUILDERCOND_H

// Athena
// Athena
#include "AthenaBaseComps/AthAlgTool.h"
#include "CxxUtils/checker_macros.h"
#include "HGTD_ReadoutGeometry/HGTD_DetectorElementCollection.h"

// Amg
#include "GeoPrimitives/GeoPrimitives.h"
// Trk
#include "TrkDetDescrInterfaces/ILayerBuilderCond.h"
#include "TrkDetDescrUtils/SharedObject.h"
#include "TrkDetDescrUtils/BinnedArray.h"
// STL
#include <vector>

class HGTD_ID;
class HGTD_DetectorManager;

namespace Trk {
  class Surface;
  class CylinderLayer;
  class DiscLayer;
  class PlaneLayer;
  class BinnedLayerMaterial;
  typedef std::pair< SharedObject<const Surface>, Amg::Vector3D > SurfaceOrderPosition;
}

/** @class HGTD_LayerBuilderCond

   The HGTD_LayerBuilderCond parses the senstive detector elments and orders them onto a
   Disc surface; no cylindrical layers are expected. 
   This implementation is based on what done in the SiLayerBuilderCond, adapted to the HGTD use case.

   */

class ATLAS_NOT_THREAD_SAFE HGTD_LayerBuilderCond :
public AthAlgTool, virtual public Trk::ILayerBuilderCond {

  public:

    /** AlgTool style constructor */
    HGTD_LayerBuilderCond(const std::string&,const std::string&,const IInterface*);

    /** Destructor */
    virtual ~HGTD_LayerBuilderCond();

    /** AlgTool initialize method */
    virtual StatusCode initialize() override;
    /** AlgTool finalize method */
    virtual StatusCode finalize() override;

    /** LayerBuilder interface method - returning Barrel-like layers */
    virtual std::unique_ptr<const std::vector<Trk::CylinderLayer*> >
    cylindricalLayers(const EventContext& ctx,
                      SG::WriteCondHandle<Trk::TrackingGeometry>& whandle) const override final;

    /** LayerBuilder interface method - returning Endcap-like layers */
    virtual std::unique_ptr<const std::vector<Trk::DiscLayer*> >
    discLayers(const EventContext& ctx,
               SG::WriteCondHandle<Trk::TrackingGeometry>& whandle) const override final;

    /** LayerBuilder interface method - returning Planar-like layers */
    virtual std::unique_ptr<const std::vector<Trk::PlaneLayer*> >
    planarLayers(const EventContext& ctx,
                 SG::WriteCondHandle<Trk::TrackingGeometry>& whandle) const override final;

    /** Name identification */
    virtual const std::string& identification() const override final;

  private:
    SG::ReadCondHandle<InDetDD::HGTD_DetectorElementCollection> retrieveHGTDdetElements(const EventContext& ctx) const;
    //!< helper method to construct HGTD materia 
    const Trk::BinnedLayerMaterial discLayerMaterial(double rMin, double rMax) const;

    //!< layer association
    void registerSurfacesToLayer( Trk::BinnedArraySpan<Trk::Surface const * const>& surfaces,const Trk::Layer& layer) const;

    static void evaluateBestBinning(std::vector<Trk::SurfaceOrderPosition>& surfaces,
                             std::vector<float>& rBins,
                             float& maxRadius,
                             std::vector<std::vector<float>>& phiBins) ;

    const HGTD_DetectorManager*           m_hgtdMgr;                        //!< the HGTD Detector Manager
    const HGTD_ID*                        m_hgtdHelper;                     //!< HGTD Id Helper
                                          
    bool                                  m_setLayerAssociation;            //!< Set Layer Association
                                          
    std::string                           m_identification;                 //!< string identification        
    
    int                                   m_rBins;                          //!< set the number of bins
    int                                   m_phiBins;                        //!< set the number of bins
    
    float                                 m_discEnvelopeR;                  //!< set disc envelope
    float                                 m_discThickness;                  //!< set disc thickness
    
    bool                                  m_runGeometryValidation;          //!< run geometry validation

    SG::ReadCondHandleKey<InDetDD::HGTD_DetectorElementCollection>
      m_HGTD_ReadKey{
        this,
        "HGTD_ReadKey",
        "HGTD_DetectorElementCollection",
        "Key of output HGTD_DetectorElementCollection for HGTD"
      };
};

inline std::unique_ptr<const std::vector<Trk::CylinderLayer*> >
HGTD_LayerBuilderCond::cylindricalLayers(const EventContext&,
                                         SG::WriteCondHandle<Trk::TrackingGeometry>& /*whandle*/) const
{
  return nullptr;
}

inline std::unique_ptr<const std::vector<Trk::PlaneLayer*> >
HGTD_LayerBuilderCond::planarLayers(const EventContext&,
                                    SG::WriteCondHandle<Trk::TrackingGeometry>& /*whandle*/) const
{
  return nullptr;
}

inline const std::string& HGTD_LayerBuilderCond::identification() const
{ return m_identification; }

#endif // HGTD_TRACKINGGEOMETRY_HGTDLAYERBUILDERCOND_H
