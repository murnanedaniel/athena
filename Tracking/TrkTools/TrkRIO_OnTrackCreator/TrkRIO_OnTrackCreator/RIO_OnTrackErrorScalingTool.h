/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

// RIO_OnTrackErrorScalingTool.h
//   Header file for class RIO_OnTrackErrorScalingTool
///////////////////////////////////////////////////////////////////
// (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////
// Wolfgang.Liebig@cern.ch
///////////////////////////////////////////////////////////////////

#ifndef TRKTOOLS_RIO_ONTRACKERRORSCALINGTOOL_H
#define TRKTOOLS_RIO_ONTRACKERRORSCALINGTOOL_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/MsgStream.h"
#include "AthenaKernel/IIOVSvc.h"

#include "TrkToolInterfaces/IRIO_OnTrackErrorScalingTool.h"
#include "EventPrimitives/EventPrimitives.h"

namespace Trk {

  /** @class RIO_OnTrackErrorScalingTool
      @brief adapts measurement errors to better reflect insufficiently
             calibrated/aligned data, using parameters from COOL.

      The RIO_OnTrackErrorScalingTool applies a linear function,

      @f$ \mathrm{error}' = a^2\cdot\mathrm{error}^2 + c^2 @f$

      to inflate measurement errors for track fitting in the reconstruction.
      It provides a handle for call-back by COOL after detector conditions
      are known. The parameters @f$a,c@f$ are retrieved from two folders
      in COOL: /Indet/TrkErrorScaling and /MUON/TrkErrorScaling.

      For performance reasons (i.e. avoid repeated re-creation of
      Trk::Covariance objects in memory) the chosen design is a
      covariance matrix factory and therefore expects
      "good cooperation" from the detector-specific Trk::IRIO_OnTrackCreator
      implementations: they need to call the associated scaling method
      from this tool at the point when covariance matrix objects
      are made.

      Note that this tool can not replace a final, precise detector
      calibration. The parameters are therefore sub-divided into
      barrel and end-cap technologies, and according to the number
      of local coordinate measurements provided by each detector,
      but currently do not allow a more finely granulated scaling.

      @author Wolfgang Liebig <http://consult.cern.ch/xwho/people/54608>
  */
  class RIO_OnTrackErrorScalingTool : public AthAlgTool,
        virtual public IRIO_OnTrackErrorScalingTool {

  public:
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////
    //! standard AlgTool constructor
    RIO_OnTrackErrorScalingTool(const std::string&,const std::string&,
		       const IInterface*);
    virtual ~RIO_OnTrackErrorScalingTool();
    //! standard AlgTool initialisation
    virtual StatusCode initialize();
    //! standard AlgTool termination
    virtual StatusCode finalize();
    //! methods indicating if the client (concrete sub-det tools) should call the scaling tool at each meas't.
    virtual bool needToScalePixel() const;
    virtual bool needToScaleSct() const; 
    virtual bool needToScaleTrt() const; 
    virtual bool needToScaleMdt() const; 
    virtual bool needToScaleTgc() const; 
    virtual bool needToScaleRpc() const; 
    virtual bool needToScaleCsc() const; 
    //! cov matrix creation from cov input and Pixel error scaling parameters
    virtual Amg::MatrixX* createScaledPixelCovariance
      (const Amg::MatrixX&, bool) const;
    //! cov matrix creation from cov input and SCT error scaling parameters 
    virtual Amg::MatrixX* createScaledSctCovariance
      (const Amg::MatrixX&, bool, double ) const;
    //! cov matrix creation from cov input and TRT error scaling parameters 
    virtual Amg::MatrixX* createScaledTrtCovariance
      (const Amg::MatrixX&, bool) const;
    //! cov matrix creation from cov input and MDT error scaling parameters 
    virtual Amg::MatrixX* createScaledMdtCovariance
      (const Amg::MatrixX&, bool) const;
    //! cov matrix creation from cov input and TGC error scaling parameters 
    virtual Amg::MatrixX* createScaledTgcCovariance
      (const Amg::MatrixX&, const ParamDefs) const;
    //! cov matrix creation from cov input and RPC error scaling parameters 
    virtual Amg::MatrixX* createScaledRpcCovariance
      (const Amg::MatrixX&, const ParamDefs) const;
    //! cov matrix creation from cov input and CSC error scaling parameters 
    virtual Amg::MatrixX* createScaledCscCovariance
      (const Amg::MatrixX&, const ParamDefs) const;

    //! dumps the current scaling parameters to MsgStream output
    virtual MsgStream& dump( MsgStream& ) const;
  private:
  ///////////////////////////////////////////////////////////////////
  // Private data members:
  ///////////////////////////////////////////////////////////////////
    bool                  m_do_pix;
    bool                  m_do_sct;
    bool                  m_do_trt;
    bool                  m_do_mdt;
    bool                  m_do_tgc;
    bool                  m_do_rpc;
    bool                  m_do_csc;
    std::vector<double>   m_scaling_pixPhi_barrel;
    std::vector<double>   m_scaling_pixPhi_endcap;
    std::vector<double>   m_scaling_pixEta_barrel;
    std::vector<double>   m_scaling_pixEta_endcap;
    std::vector<double>   m_scaling_sct_barrel;
    std::vector<double>   m_scaling_sct_endcap;
    std::vector<double>   m_scaling_trt_barrel;
    std::vector<double>   m_scaling_trt_endcap;
    std::vector<double>   m_scaling_mdt_barrel;
    std::vector<double>   m_scaling_mdt_endcap;
    std::vector<double>   m_scaling_tgcPhi;
    std::vector<double>   m_scaling_tgcEta;
    std::vector<double>   m_scaling_rpcPhi;
    std::vector<double>   m_scaling_rpcEta;
    std::vector<double>   m_scaling_cscPhi;
    std::vector<double>   m_scaling_cscEta;

    //! conditions data handling: COOL folder name for InDet scaling paremeters
    const std::string m_idFolder;
    //! conditions data handling: COOL folder name for Muon scaling paremeters
    const std::string m_muonFolder;
    //! map holding the scaling parameters for each detector region
    typedef std::map<std::string,std::vector<double>*> ParamMap;
    /** @brief map holding the scaling parameters for each detector region

        The parammap is used at call-back/initialisation to avoid map
        search at every scaling call. It maps directly to the
        references of the local vector variables!
     */
    ParamMap m_parammap;
    //! internal code structuring: formatted output for scaling parameters
    const std::string makeInfoString(const std::string&,const bool,
                                     const std::vector<double>&) const;
    //! conditions data handling: find associated parameters for a given detector
    void registerParameters(bool&,
                            std::vector<double>*,
                            const std::string&);

    //! internal code structuring: error scaling for 2-dim measurements
    void scale2by2(Amg::MatrixX&, 
                   const std::vector<double>&,
                   const std::vector<double>&) const;

    //! conditions data handling: call-back entry to re-set scaling parameters when an IOV is new or changed.
    StatusCode callback( IOVSVC_CALLBACK_ARGS );
  };

  ///////////////////////////////////////////////////////////////////
  // Inline methods:
  ///////////////////////////////////////////////////////////////////
  inline bool RIO_OnTrackErrorScalingTool::needToScalePixel() const
  {
    return m_do_pix;
  }
  inline bool RIO_OnTrackErrorScalingTool::needToScaleSct() const
  {
    return m_do_sct;
  }
  inline bool RIO_OnTrackErrorScalingTool::needToScaleTrt() const
  {
    return m_do_trt;
  }
  inline bool RIO_OnTrackErrorScalingTool::needToScaleMdt() const
  {
    return m_do_mdt;
  }
  inline  bool RIO_OnTrackErrorScalingTool::needToScaleTgc() const
  {
    return m_do_tgc;
  }
  inline bool RIO_OnTrackErrorScalingTool::needToScaleRpc() const
  {
    return m_do_rpc;
  }
  inline bool RIO_OnTrackErrorScalingTool::needToScaleCsc() const
  {
    return m_do_csc;
  }

  MsgStream& operator << (MsgStream&, const RIO_OnTrackErrorScalingTool&);
  ///////////////////////////////////////////////////////////////////
  // Overload of << operator MsgStream
  ///////////////////////////////////////////////////////////////////
  inline MsgStream& operator    << 
    (MsgStream& sl,const RIO_OnTrackErrorScalingTool& se)
    { 
      return se.dump(sl); 
    }

} // end of namespace

#endif
