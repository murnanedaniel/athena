/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef CALOCLUSTERCORRECTION_CALOSWDEADOTX_PS_H
#define CALOCLUSTERCORRECTION_CALOSWDEADOTX_PS_H
/********************************************************************

NAME:     CaloSwDeadOTX_ps.h
PACKAGE:  offline/Calorimeter/CaloClusterCorrection
          
********************************************************************/

#include "GaudiKernel/ToolHandle.h"
#include "CaloInterface/ICaloAffectedTool.h"
#include "CaloClusterCorrection/CaloClusterCorrectionCommon.h"
#include "CaloConditions/Array.h"

class CaloSwDeadOTX_ps : public CaloClusterCorrectionCommon
{
public:
  // constructor 
  CaloSwDeadOTX_ps (const std::string& type,
                     const std::string& name,
                     const IInterface* parent);

  /**
   * @brief Virtual function for the correction-specific code.
   * @param cluster The cluster to correct.
   *                It is updated in place.
   * @param elt     The detector description element corresponding
   *                to the cluster location.
   * @param eta     The @f$\eta@f$ coordinate of the cluster, in this sampling.
   * @param adj_eta The @f$\eta@f$ adjusted for
   *                any shift between the actual and nominal coordinates.
   *                (This is shifted back to the nominal coordinate system.)
   * @param phi     The @f$\phi@f$ coordinate of the cluster, in this sampling.
   * @param adj_phi The @f$\phi@f$ adjusted for
   *                any shift between the actual and nominal coordinates.
   *                (This is shifted back to the nominal coordinate system.)
   * @param samp    The calorimeter sampling we're examining.
   *                This is a sampling code as defined by
   *                @c CaloSampling::CaloSample; i.e., it has both
   *                the calorimeter region and sampling encoded.
   */
  virtual void makeTheCorrection (xAOD::CaloCluster* cluster,
                                  const CaloDetDescrElement* elt,
                                  float eta,
                                  float adj_eta,
                                  float phi,
                                  float adj_phi,
                                  CaloSampling::CaloSample samp) const;

  virtual StatusCode initialize();

 private:

  //  std::string m_affectedTool;
  ToolHandle<ICaloAffectedTool> m_affectedTool;
  
  CaloRec::Array<3>     m_correction;
  CaloRec::Array<2>     m_sampling_depth;
  float m_eta_start_crack;
  float m_eta_end_crack;
  float m_etamax;
  bool m_use_raw_eta;
  

};

#endif

