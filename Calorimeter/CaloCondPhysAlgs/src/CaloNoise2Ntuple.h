/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

// CaloNoise2Ntuple.h
//

#ifndef CALOCONDPHYSALGS_CALONOISE2NTUPLE_H
#define CALOCONDPHYSALGS_CALONOISE2NTUPLE_H

#include <string>

// Gaudi includes

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "CaloIdentifier/CaloIdManager.h"
#include "CaloDetDescr/CaloDetDescrManager.h"
#include "CaloIdentifier/CaloCell_ID.h"

#include "GaudiKernel/ITHistSvc.h"
#include "TTree.h"

#include "StoreGate/ReadCondHandleKey.h"

class CaloNoise;

class CaloNoise2Ntuple : public AthAlgorithm {

  public:
    //Gaudi style constructor and execution methods
    /** Standard Athena-Algorithm Constructor */
    CaloNoise2Ntuple(const std::string& name, ISvcLocator* pSvcLocator);
    /** Default Destructor */
    virtual ~CaloNoise2Ntuple();
    
    /** standard Athena-Algorithm method */
    virtual StatusCode          initialize() override;
    /** standard Athena-Algorithm method */
    virtual StatusCode          execute() override;
    /** standard Athena-Algorithm method */
    virtual StatusCode          finalize() override;
    /** standard Athena-Algorithm method */
    virtual StatusCode          stop() override;
    
  private:

  //---------------------------------------------------
  // Member variables
  //---------------------------------------------------
  ITHistSvc* m_thistSvc;

  const CaloCell_ID*       m_calo_id;

  SG::ReadCondHandleKey<CaloNoise> m_totalNoiseKey
    { this, "TotalNoiseKey", "totalNoise", "SG key for total noise" };
  SG::ReadCondHandleKey<CaloNoise> m_elecNoiseKey
    { this, "ElecNoiseKey", "electronicNoise", "SG key for electronic noise" };
  SG::ReadCondHandleKey<CaloNoise> m_pileupNoiseKey
    { this, "PileupNoiseKey", "pileupNoise", "SG key for pileup noise" };

  std::string m_treeName;

  int m_iCool;
  int m_SubHash;
  int m_Hash;
  int m_OffId;
  float m_eta;
  float m_phi;
  int m_layer;
  int m_Gain;
  float m_noise;
  float m_elecNoise;
  float m_pileupNoise; 
  TTree* m_tree;

  int m_runNumber;
  int m_lumiBlock;

};
#endif
