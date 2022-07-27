/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EGAMMATOOLS_EGAMMAOQFLAGSBUILDER_H
#define EGAMMATOOLS_EGAMMAOQFLAGSBUILDER_H
/**
  @class egammaOQFlagsBuilder
  egamma Object Quality flags data object builder :
    - This tool checks if any cell of the cluster associated to the egamma
object is affected by a detector problem: non nominal or dead high voltage,
readout problems, missing FEBs, high quality factor, timing, etc.... If this is
the case, then a bit corresponding to a specific problem is filled ( see
egammaEvent/egammaEvent/egammaPIDdefs.h for bits definition). Most of the
informations are given separately for each layer of the EM calorimeter. They are
also separately stored for the core (3x3 central cells in the middle layer of
the EM calorimeter) and the edge (all the other cells) of the cluster.
  @author Frederic Derue derue@lpnhe.in2p3.fr
  @author Francesco Polci polci@lpsc.in2p3.fr
  @author Quentin Buat quentin.buat@lpsc.in2p3.fr
*/

#include "egammaInterfaces/IegammaOQFlagsBuilder.h"

#include "AthenaBaseComps/AthAlgTool.h"

#include "CaloInterface/ICaloAffectedTool.h"
#include "CaloIdentifier/CaloCell_ID.h"
#include "CaloIdentifier/LArEM_ID.h"
#include "CaloUtils/CaloCellList.h"
#include "GaudiKernel/EventContext.h"
#include "GaudiKernel/ToolHandle.h"

#include "Identifier/HWIdentifier.h"
#include "LArCabling/LArOnOffIdMapping.h"
#include "LArRecConditions/LArBadChannelCont.h"
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/ReadHandleKey.h"
#include "CaloConditions/CaloAffectedRegionInfoVec.h"

#include "xAODCaloEvent/CaloClusterFwd.h"
#include "xAODEgamma/EgammaFwd.h"

class CaloCellContainer;
class HWIdentifier;
class LArEM_ID;
class CaloCell_ID;

class egammaOQFlagsBuilder final
  : public AthAlgTool
  , virtual public IegammaOQFlagsBuilder
{
public:
  /** @brief Default constructor*/
  egammaOQFlagsBuilder(const std::string& type,
                       const std::string& name,
                       const IInterface* parent);

  /** @brief Destructor*/
  ~egammaOQFlagsBuilder();
  /** @brief initialize method*/
  StatusCode initialize();
  /** @brief standard execute method */
  virtual StatusCode execute(const EventContext& ctx,
                             xAOD::Egamma& egamma) const;
  /** @brief finalize method*/
  StatusCode finalize();

private:
  /** Handle to bad-channel CDO */
  SG::ReadCondHandleKey<LArBadChannelCont> m_bcContKey{
    this,
    "LArBadChannelKey",
    "LArBadChannel",
    "Key of the LArBadChannelCont CDO"
  };

  SG::ReadCondHandleKey<CaloAffectedRegionInfoVec> m_affKey{
    this,
    "LArAffectedRegionKey",
    "LArAffectedRegionInfo",
    "SG key for affected regions cond object"
  };

  ToolHandle<ICaloAffectedTool> m_affectedTool{ this,
                                                "affectedTool",
                                                "CaloAffectedTool",
                                                "CaloAffectedTool" };

  const LArEM_ID* m_emHelper;
  const CaloCell_ID* m_calocellId;

  SG::ReadHandleKey<CaloCellContainer> m_cellsKey{
    this,
    "CellsName",
    "AllCalo",
    "Names of container which contain cells"
  };

  Gaudi::Property<double> m_QCellCut{ this, "QCellCut", 4000. };
  Gaudi::Property<double> m_QCellHECCut{ this, "QCellHECCut", 60000. };
  Gaudi::Property<double> m_QCellSporCut{ this, "QCellSporCut", 4000. };
  Gaudi::Property<double> m_LArQCut{ this, "LArQCut", 0.8 };
  Gaudi::Property<double> m_TCut{ this, "TCut", 10.0 };
  Gaudi::Property<double> m_TCutVsE{ this, "TCutVsE", 2.0 };
  Gaudi::Property<double> m_RcellCut{ this, "RcellCut", 0.8 };
};

#endif

