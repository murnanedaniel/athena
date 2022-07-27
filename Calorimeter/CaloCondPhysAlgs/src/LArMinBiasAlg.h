/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/**
 *  @file TheLArMinBiasAlg.h
 */

#ifndef CALOCONDPHYSALGS_LARMINBIASALG_H
#define CALOCONDPHYSALGS_LARMINBIASALG_H

#include <string>

// Gaudi includes

#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "StoreGate/ReadHandleKey.h"
#include "xAODEventInfo/EventInfo.h"

#include "CaloDetDescr/CaloDetDescrManager.h"


#include "LArCabling/LArOnOffIdMapping.h"
#include "LArRawConditions/LArMCSym.h"


class TTree;
class LArEM_ID;
class CaloCell_ID;
class ITHistSvc;


#define MAX_SYM_CELLS 3000

  class LArMinBiasAlg : public AthAlgorithm {
  public:
    //Gaudi style constructor and execution methods
    /** Standard Athena-Algorithm Constructor */
    LArMinBiasAlg(const std::string& name, ISvcLocator* pSvcLocator);
    /** Default Destructor */
    ~LArMinBiasAlg();
    
    virtual StatusCode  initialize() override;
    virtual StatusCode  execute() override;
    virtual StatusCode  finalize() override;
    virtual StatusCode  stop() override;

    
  private:

     void        fillNtuple();
     void        addCell(int index, double e1, double e2, double wt=1. );

  //---------------------------------------------------
  // Member variables
  //---------------------------------------------------
  int m_datasetID_lowPt;
  int m_datasetID_highPt;
  double m_weight_lowPt;
  double m_weight_highPt;
  SG::ReadCondHandleKey<LArMCSym> m_mcSymKey
  { this, "MCSymKey", "LArMCSym", "SG Key of LArMCSym object" };
  SG::ReadCondHandleKey<LArOnOffIdMapping> m_cablingKey{this,"CablingKey","LArOnOffIdMap","SG Key of LArOnOffIdMapping object"};

  SG::ReadCondHandleKey<CaloDetDescrManager> m_caloMgrKey { this
      , "CaloDetDescrManager"
      , "CaloDetDescrManager"
      , "SG Key for CaloDetDescrManager in the Condition Store" };

  const LArEM_ID*        m_larem_id = nullptr;
  const CaloCell_ID*       m_calo_id = nullptr;
  std::vector<double> m_eCell;
  

  ITHistSvc* m_thistSvc = nullptr;
  TTree* m_tree = nullptr;
  int m_nevt_total = 0;
  int m_n1 = 0;
  int m_n2 = 0;

// FIXME   Total maximum array size for ntuple hardcoded... not very nice
   int m_nsymcell = 0;
   double m_nevt[MAX_SYM_CELLS]{};
   int m_layer[MAX_SYM_CELLS]{};
   int m_region[MAX_SYM_CELLS]{};
   int m_identifier[MAX_SYM_CELLS]{};
   int m_ieta[MAX_SYM_CELLS]{};
   float m_eta[MAX_SYM_CELLS]{};
   float m_phi[MAX_SYM_CELLS]{};
   double m_average[MAX_SYM_CELLS]{};
   double m_rms[MAX_SYM_CELLS]{};
   double m_offset[MAX_SYM_CELLS]{};


  struct CellInfo {
      int layer;
      int region;
      int ieta;
      float eta;
      float phi;
      Identifier identifier;
      double nevt;
      double average;
      double rms;
      double offset;
  };
  std::vector<CellInfo> m_CellList;
  std::vector<int> m_symCellIndex;
  float m_first;
  int m_ncell = 0;
  SG::ReadHandleKey<xAOD::EventInfo> m_eventInfoKey{this,"EvtInfo", "EventInfo", "EventInfo name"};
  };
#endif
