/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#ifndef HIEVENTUTILS_HITOWERWEIGHTTOOL_H
#define HIEVENTUTILS_HITOWERWEIGHTTOOL_H

#include "HIEventUtils/IHITowerWeightTool.h"
#include "AthenaBaseComps/AthAlgTool.h"

#include <TH3F.h>
#include <TFile.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>


////////////////////////////////////////////////////////////
// Tool that gets the tower weights for by tower by tower
// difference in response during the HI jet reconstruction.
////////////////////////////////////////////////////////////

class HITowerWeightTool : public extends<AthAlgTool, IHITowerWeightTool>
{
  public:
    HITowerWeightTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~HITowerWeightTool() = default;
    virtual StatusCode initialize() override;

    virtual float getEtaPhiResponse(float eta, float phi, const EventContext& ctx) const override;
    virtual float getEtaPhiOffset(float eta, float phi, const EventContext& ctx) const override;
    virtual float getWeight(float eta, float phi, int sampling) const override;
    virtual float getWeightEta(float eta, float phi, int sampling) const override;
    virtual float getWeightPhi(float eta, float phi, int sampling) const override;
    virtual float getWeightMag(float eta, float phi, int sampling) const override;
    
  private:
    Gaudi::Property<bool> m_applycorrection{this, "ApplyCorrection", true , "If false, unit weigts are applied"};
    Gaudi::Property<std::string> m_inputFile{this, "InputFile", "cluster.geo.HIJING_2018.root","File containing cluster geometric moments."};
    Gaudi::Property<std::string> m_configDir{this, "ConfigDir", "HIJetCorrection/","Directory containing configuration file."};

    TH3F* m_h3W;
    TH3F* m_h3Eta;
    TH3F* m_h3Phi;
    TH3F* m_h3Mag;
    TH3F* m_h3EtaPhiResponse;
    TH3F* m_h3EtaPhiOffset;
    std::map<unsigned int, int> m_runMap;

    int getRunIndex(const EventContext& ctx) const;

};

#endif
