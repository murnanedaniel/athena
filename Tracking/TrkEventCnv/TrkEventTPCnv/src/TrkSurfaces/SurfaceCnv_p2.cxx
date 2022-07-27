/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//-----------------------------------------------------------------------------
//
// file:   SurfaceCnv_p2.cxx
//
// ----------------------------------------------------------------------------

#include "TrkSurfaces/Surface.h"
#include "TrkSurfaces/StraightLineSurface.h"
#include "TrkSurfaces/DiscSurface.h"
#include "TrkSurfaces/PlaneSurface.h"
#include "TrkSurfaces/CylinderSurface.h"
#include "TrkSurfaces/PerigeeSurface.h"
#include "TrkDistortedSurfaces/SaggedLineSurface.h"
#include "TrkEventTPCnv/TrkSurfaces/SurfaceCnv_p2.h"
#include "TrkEventTPCnv/helpers/EigenHelpers.h"
#include "CxxUtils/checker_macros.h"

template <class SURFACE>
SURFACE* SurfaceCnv_p2<SURFACE>::createTransient( const Trk::Surface_p2 * persObj,MsgStream& ){
  SURFACE* surface=nullptr;
  if (persObj->m_transform.empty()) {
    // det element surface
    Identifier id =  Identifier32(persObj->m_associatedDetElementId);
    const SURFACE* detSurf =  static_cast<const SURFACE*>(m_eventCnvTool->getSurface(id));
    // In this case, we return a surface object held by detdescr.
    // We need to cast away const to match the createTransient interface.
    // This is ugly, but in practice should be ok.
    // FIXME: can this be done better?
    SURFACE* surface_nc ATLAS_THREAD_SAFE = const_cast<SURFACE*>(detSurf);
    surface= surface_nc; // Needed to fulfill interface...
  } else {
    // Not Det element surface, so need to create surface & fill transform
    Amg::Transform3D transform;
    EigenHelpers::vectorToEigenTransform3D( persObj->m_transform, transform );
    surface=new SURFACE(transform);
  }  
  return surface;
}

template <class SURFACE>
void SurfaceCnv_p2<SURFACE>::persToTrans( const Trk::Surface_p2 * , SURFACE * , MsgStream &) {
  throw std::runtime_error("SurfaceCnv_p2::persToTrans shouldn't be called any more!");
  // The reason for this is that persToTrans relies on an object being created already. 
  // But for det element surfaces, we will return one which already exists and belongs to the detector elements.
}

template <class	SURFACE>
void SurfaceCnv_p2<SURFACE>::transToPers( const SURFACE         * transObj,
                                       Trk::Surface_p2 * persObj,
                                       MsgStream & log)
{
  if (not transObj){
    log<<MSG::WARNING<<"SurfaceCnv_p2<SURFACE>::transToPers - null pointer for transient object."<<endmsg;
    return;
  }
  persObj->m_associatedDetElementId = transObj->associatedDetectorElementIdentifier().get_identifier32().get_compact();
  if (transObj->isFree() ) { 
    // if this is a free surface, write it out 'as is' 
    // - otherwise we're relying on the fact that the Identifier is there (and valid!) to reconstruct it on reading
    EigenHelpers::eigenTransform3DToVector( *transObj->cachedTransform(),persObj->m_transform );
    
  } else {
    if (!transObj->associatedDetectorElementIdentifier().is_valid()) {
       log<<MSG::WARNING<<"SurfaceCnv_p2<SURFACE>::transToPers - invalid detector element for non-free surface: "<<*transObj<<endmsg;
    }
  }
  persObj->m_surfaceType = static_cast<uint8_t>(transObj->type());
}

// Specialisation for SaggedLineSurface - hopefully this will only ever be a DetEl surface, otherwise we need a new persistent type.
template <>
void SurfaceCnv_p2<Trk::SaggedLineSurface>::transToPers( const Trk::SaggedLineSurface         * transObj,
                                       Trk::Surface_p2 * persObj,
                                       MsgStream &log)
{
  persObj->m_associatedDetElementId = transObj->associatedDetectorElementIdentifier().get_identifier32().get_compact();
  persObj->m_surfaceType = static_cast<uint8_t>(transObj->type());
  if (transObj->isFree() ) {
    log << MSG::WARNING<<"Not expecting these to be Free!"<<endmsg;
  }   
}

template class SurfaceCnv_p2<Trk::ConeSurface>;
template class SurfaceCnv_p2<Trk::CylinderSurface>;
template class SurfaceCnv_p2<Trk::DiscSurface>; 
template class SurfaceCnv_p2<Trk::PerigeeSurface>;
template class SurfaceCnv_p2<Trk::PlaneSurface>;
template class SurfaceCnv_p2<Trk::StraightLineSurface>;
template class SurfaceCnv_p2<Trk::SaggedLineSurface>;
