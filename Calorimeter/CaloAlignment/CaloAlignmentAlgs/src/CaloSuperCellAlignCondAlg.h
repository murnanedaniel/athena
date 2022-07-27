/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef CALOALIGNMENTALGS_CALOSUPERCELLALIGNCONDALG_H
#define CALOALIGNMENTALGS_CALOSUPERCELLALIGNCONDALG_H

#include "AthenaBaseComps/AthAlgorithm.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/WriteCondHandleKey.h"
#include "GaudiKernel/ToolHandle.h"

#include "CaloDetDescr/CaloDetDescrManager.h"
#include "CaloDetDescr/ICaloSuperCellIDTool.h"

/**
 * @class CaloSuperCellAlignCondAlg
 *
 * @brief Condition Algorithm for making CaloSuperCellDetDescrManager condition object
 *
 **/

class CaloSuperCellAlignCondAlg final : public AthAlgorithm
{
 public:
  using AthAlgorithm::AthAlgorithm;
  virtual ~CaloSuperCellAlignCondAlg() = default;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

 private:
  SG::ReadCondHandleKey<CaloDetDescrManager>  m_readCaloMgrKey {this
      , "CaloDetDescrManager"
      , "CaloDetDescrManager"
      , "SG key of the resulting CaloDetDescrManager" };

  SG::WriteCondHandleKey<CaloSuperCellDetDescrManager>  m_writeCaloSuperCellMgrKey {this
      , "CaloSuperCellDetDescrManager"
      , "CaloSuperCellDetDescrManager"
      , "SG key of the resulting CaloSuperCellDetDescrManager" };

  ToolHandle<ICaloSuperCellIDTool> m_scidTool {this
      , "CaloSuperCellIDTool"
      , "CaloSuperCellIDTool"
      , "Offline / SuperCell ID mapping tool" };
};

#endif
