/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


// ********************************************************************
//
// NAME:     L1CaloMonitoringCaloTool.h
// PACKAGE:  TrigT1CaloCalibTools
//
// AUTHOR:   Peter Faulkner
//
// ********************************************************************
#ifndef TRIGT1CALOCALIBTOOLS_L1CALOMONITORINGCALOTOOL_H
#define TRIGT1CALOCALIBTOOLS_L1CALOMONITORINGCALOTOOL_H

#include <string>
#include <vector>

#include "GaudiKernel/ToolHandle.h"

#include "AthenaBaseComps/AthAlgTool.h"
#include "TrigT1CaloCalibToolInterfaces/IL1CaloMonitoringCaloTool.h"

class IInterface;
class StatusCode;

class Identifier;
class CaloLVL1_ID;

namespace LVL1 {

  class IL1CaloCells2TriggerTowers;

  class L1CaloMonitoringCaloTool: public IL1CaloMonitoringCaloTool, public AthAlgTool
  {

   public:
  
    L1CaloMonitoringCaloTool(const std::string & type, const std::string & name,
		             const IInterface* parent);
    

    virtual ~L1CaloMonitoringCaloTool();

    virtual StatusCode initialize();
    virtual StatusCode finalize();

    StatusCode loadCaloCells();
    float et(const Identifier& ttid) const;
    float caloQuality(const Identifier& ttid) const;

   private:

    int towerIndex(const Identifier& ttId) const;
    int region(int index) const;
    int etaBin(int index) const;

    ToolHandle<LVL1::IL1CaloCells2TriggerTowers> m_cells2tt;
    const CaloLVL1_ID* m_lvl1Helper;

    std::string m_caloCellContainerName;

    std::vector<float> m_energySums;
    std::vector<float> m_quality;
    std::vector<float> m_denom;
    std::vector<float> m_sinTh;
    std::vector<unsigned int> m_cellIds;
    std::vector<int> m_ttIdx;
    unsigned int m_maxCells;

    int m_events;
    int m_lastRun;
    int m_lastEvent;
    int m_sideOffset;
    int m_layerOffset;
    std::vector<int> m_binOffset;
    std::vector<int> m_indexOffset;
    std::vector<int> m_etaShift;
  
    static const int s_maxTowers = 7168;
    static const int s_nregions  = 4;
    static const int s_nsinThBins = 33;
  };

} // end namespace

#endif
