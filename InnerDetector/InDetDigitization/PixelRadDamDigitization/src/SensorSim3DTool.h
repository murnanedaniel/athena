/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/

///////////////////////////////////////////////////////////////////
// SensorSim3DTool.h
//   Header file for class SensorSim3DTool
///////////////////////////////////////////////////////////////////
// (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////
//
// Configurable Parameters
//
// -numberOfSteps	Geant4:number of steps for PixelPlanar //ALEX: Pixel3D?
// -numberOfCharges	Geant4:number of charges for PixelPlanar //ALEX: Pixel3D?
//
//////////////////////////////////////////////////////////////////

#ifndef PIXELDIGITIZATION_SensorSim3DTool_H
#define PIXELDIGITIZATION_SensorSim3DTool_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "StoreGate/WriteHandle.h"

#include "SensorSimTool.h"
#include "RadDamageUtil.h"

namespace RadDam{

class SensorSim3DTool : public SensorSimTool {

  public:
    SensorSim3DTool( const std::string& type, const std::string& name,const IInterface* parent);
    virtual StatusCode initialize();
    virtual StatusCode finalize();
    virtual ~SensorSim3DTool();

    virtual StatusCode induceCharge(const TimedHitPtr<SiHit> &phit, SiChargedDiodeCollection& chargedDiodes, const InDetDD::SiDetectorElement &Module, const InDetDD::PixelModuleDesign &p_design, std::vector< std::pair<double,double> > &trfHitRecord, std::vector<double> &initialConditions);  
    //Apply slim edge inefficiencies for IBL sensors
    virtual StatusCode applySlimEdges( double &energyPerStep, double &eta_drifted);

    //Maps
    std::map<std::pair<int, int>, TH3F* > ramoPotentialMap;
    std::map<std::pair<int, int>, TH1F*> eFieldMap;
    std::map<std::pair<int, int>, TH2F*> distanceMap_e;
    std::map<std::pair<int, int>, TH2F*> distanceMap_h;
    std::map<std::pair<int, int>, TH1F*> timeMap_e;
    std::map<std::pair<int, int>, TH1F*> timeMap_h;
    std::map<std::pair<int, int>, TH2F*> lorentzMap_e;
    std::map<std::pair<int, int>, TH2F*> lorentzMap_h;

    ToolHandle<RadDamageUtil>                     m_radDamageUtil;
  
  private:
    SensorSim3DTool();

    int    m_numberOfSteps;
    int    m_numberOfCharges;  
    double m_diffusionConstant;

    bool m_doRadDamage;
    double m_fluence; //eventually, this should be pulled from the conditions.

    std::vector<double> fluence_layers;
    std::map<std::pair<int, int>, double> fluence_layersMaps;
    double m_trappingTimeElectrons;
    double m_trappingTimeHoles;
};

}
#endif // PIXELDIGITIZATION_SensorSim3DTool_H
