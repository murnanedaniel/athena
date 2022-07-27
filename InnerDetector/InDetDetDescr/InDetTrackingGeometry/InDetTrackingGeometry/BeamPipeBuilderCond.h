/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// BeamPipeBuilderCond.h, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#ifndef INDETTRACKINGGEOMETRY_BEAMPIPEBUILDERCOND_H
#define INDETTRACKINGGEOMETRY_BEAMPIPEBUILDERCOND_H

// Athena
#include "AthenaBaseComps/AthAlgTool.h"
// Trk
#include "TrkDetDescrInterfaces/ILayerBuilderCond.h"
// STL
#include <vector>
#include <utility>

namespace Trk {
  class CylinderLayer;
  class DiscLayer;
  class PlaneLayer;
}

class BeamPipeDetectorManager;

namespace InDet {

  /** @class BeamPipeBuilderCond
      Simple LayerBuilder for the BeamPipe, 
      can be configured through jobOptions:
        - radius
        - halflength
        - thickness
        - MaterialProperties
        
      later on the slight shift/rotation of the BeamPipe can be implemented
      - make a binding to the database afterwards
       
      @author Andreas.Salzburger@cern.ch 
    */
  class BeamPipeBuilderCond : public AthAlgTool,
                          virtual public Trk::ILayerBuilderCond {
    
                            
    public:
      /** AlgTool style constructor */
      BeamPipeBuilderCond(const std::string&,const std::string&,const IInterface*);
      /** Destructor */
      virtual ~BeamPipeBuilderCond();

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
      virtual const std::string& identification() const override;      
          
    private:
      
      bool                                          m_beamPipeFromDb;    //!< steer beam pipe parameters from DataBase
      const BeamPipeDetectorManager*                m_beamPipeMgr;       //!< the beam pipe manager
      std::string                                   m_beamPipeMgrName;   //!< the name of the beam pipe manager to be configured
      double                                        m_beamPipeEnvelope;  //!< radial envelope when taking the Top volume radius
      
      double                                        m_beamPipeOffsetX;    //!< beam pipe offset in x
      double                                        m_beamPipeOffsetY;    //!< beam pipe offset in y        
      double                                        m_beamPipeRadius;     //!< radius of the beam pipe
      double                                        m_beamPipeHalflength; //!< halflength of the beampipe
      double                                        m_beamPipeThickness;  //!< thickness of the beam pipe
      double                                        m_beamPipeX0;         //!< X0 of the beam pipe
      double                                        m_beamPipeL0;         //!< X0 of the beam pipe
      //double                                        m_beamPipedEdX;       //!< dEdX of the beam pipe
      double                                        m_beamPipeA;          //!< averageA of the beam pipe
      double                                        m_beamPipeZ;          //!< averageZ of the beam pipe
      double                                        m_beamPipeRho;        //!< averageRho of the beam pipe

      unsigned int                                  m_beamPipeBinsZ;      //!< number of bins in the beam pipe

      std::string                                   m_identification;     //!< string identification
      
      
  };

  inline std::unique_ptr<const std::vector<Trk::DiscLayer*> >
  BeamPipeBuilderCond::discLayers(const EventContext&,
                                  SG::WriteCondHandle<Trk::TrackingGeometry>& /*whandle*/) const
  {
    return nullptr;
  }

  inline std::unique_ptr<const std::vector<Trk::PlaneLayer*> >
  BeamPipeBuilderCond::planarLayers(const EventContext&,
                                    SG::WriteCondHandle<Trk::TrackingGeometry>& /*whandle*/) const
  {
    return nullptr;
  }

 inline const std::string& BeamPipeBuilderCond::identification() const
 { return m_identification; } 
 
  
} // end of namespace


#endif // INDETTRACKINGGEOMETRY_BEAMPIPEBUILDERCOND_H


