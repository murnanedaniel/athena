/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
//  Header file for class  MMClusterOnTrackCreator
///////////////////////////////////////////////////////////////////
// (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////
// Interface for CscClusterOnTrack production
///////////////////////////////////////////////////////////////////
// Version 1.0 20/07/2004 
///////////////////////////////////////////////////////////////////

#ifndef MMClusterOnTrackCreator_H
#define MMClusterOnTrackCreator_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "MuonRecToolInterfaces/IMuonClusterOnTrackCreator.h"
#include "TrkToolInterfaces/IRIO_OnTrackErrorScalingTool.h"
#include "MuonRIO_OnTrack/MuonClusterOnTrack.h"

#include "MuonReadoutGeometry/MuonDetectorManager.h"
#include "TrkPrepRawData/PrepRawDataCLASS_DEF.h"
#include "TrkParameters/TrackParameters.h"

#include "MuonPrepRawData/CscStripPrepDataContainer.h"
#include "MuonIdHelpers/MuonIdHelperTool.h"




class ICscStripFitter;
class ICscClusterFitter;
class ICscClusterUtilTool;
class CscIdHelper;

namespace Muon {


  /** @class MMClusterOnTrackCreator
      @brief Interface for the reconstruction to calibration and alignment corrections. It should be used by 
             reconstruction and pattern recognition to create Muon::MuonClusterOnTrack objects (s).

       It offers several interfaces:
       - Create new Muon::MuonClusterOnTrack from a Trk::PrepRawData and a predicted Trk::TrackParameter.       
         @code const MuonClusterOnTrack* correct ( const Trk::PrepRawData& RIO, const Trk::TrackParameters& tp) const @endcode
       - Create new Muon::MuonClusterOnTrack from a Trk::PrepRawData and a prediction of the global position and direction.
         @code createRIO_OnTrack(const Trk::PrepRawData& ROP, const Trk::GlobalPosition& GP, const Trk::GlobalDirection GD) const @endcode
       - Create new Muon::MuonClusterOnTrack from a Trk::PrepRawData and a prediction intersect position of the muon with the 
         measurement surface.
         Kept for legacy with interface
         @code createRIO_OnTrack(const Trk::PrepRawData& RIO, const Trk::GlobalPosition& GP) const @endcode
       
       JobOptions Flags:
   */
 class MMClusterOnTrackCreator : public AthAlgTool, virtual public IMuonClusterOnTrackCreator {
  public:
    
    MMClusterOnTrackCreator(const std::string&,const std::string&,const IInterface*);
    virtual ~MMClusterOnTrackCreator()=default;
    virtual StatusCode initialize();
 
    /** @brief Create new Muon::MuonClusterOnTrack from a Trk::PrepRawData and a predicted Trk::TrackParameter. 
	@param RIO Trk::PrepRawData object to be calibrated
	@param GP  Predicted intersect position of the muon with the measurement plane 
	@return a pointer to a new Muon::MuonClusterOnTrack object, zero if calibration failed.
	The ownership of the new Muon::MuonClusterOnTrack is passed to the client calling the tool
    */       
    virtual const MuonClusterOnTrack* createRIO_OnTrack(const Trk::PrepRawData& RIO,
                                                        const Amg::Vector3D& GP) const;

    /** @brief Create new Muon::MuonClusterOnTrack from a Trk::PrepRawData and a prediction of the global position and direction.
	It is only implemented for the CSCs, for RPC and TGC Trk::PrepRawData the result is the same as for the routine without the direction. 
	@param RIO Trk::PrepRawData object to be calibrated
	@param GP  Predicted intersect position of the muon with the measurement plane 
	@param GD  Predicted direction at the intersect position of the muon with the measurement plane 
	@return a pointer to a new Muon::MuonClusterOnTrack object, zero if calibration failed.
	The ownership of the new Muon::MuonClusterOnTrack is passed to the client calling the tool
    */
    virtual const MuonClusterOnTrack* createRIO_OnTrack(const Trk::PrepRawData& RIO,
                                                        const Amg::Vector3D& GP,
                                                        const Amg::Vector3D& GD) const;
 
    /** @brief Create new Muon::MuonClusterOnTrack from a Trk::PrepRawData and the predicted Trk::TrackParameter at the measurement surface. 
	@param RIO Trk::PrepRawData object to be calibrated
	@param TP  Predicted Trk::TrackParameter at the measurement surface
	@return a pointer to a new Muon::MuonClusterOnTrack object, zero if calibration failed.
	The ownership of the new Muon::MuonClusterOnTrack is passed to the client calling the tool
    */
    virtual const MuonClusterOnTrack* correct(const Trk::PrepRawData& RIO,const Trk::TrackParameters& TP) const; 
  
  private:

    ToolHandle<MuonIdHelperTool> m_muonIdHelperTool;

 }; // end of class def
}  //  namespace Muon
#endif // MMClusterOnTrackCreator_H
