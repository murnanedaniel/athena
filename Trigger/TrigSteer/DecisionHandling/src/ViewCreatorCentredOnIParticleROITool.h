/*
Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef DECISIONHANDLING_VIEWCREATORCENTREDONIPARTICLEROITOOL_H
#define DECISIONHANDLING_VIEWCREATORCENTREDONIPARTICLEROITOOL_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "StoreGate/WriteHandleKey.h"
#include "DecisionHandling/IViewCreatorROITool.h"
#include "TrigSteeringEvent/TrigRoiDescriptorCollection.h"
#include "xAODMuon/Muon.h"

/**
 * @class ViewCreatorCentredOnIParticleROITool
 * Creates a new ROI centred on an object deriving from xAOD::IParticle.
 *
 * Stores this new ROI in the output container, and links it to the Decision Object
 *
 * The new EventView spawned by the parent EventViewCreatorAlgorithm of this tool will process in this new ROI.
 *
 * In the majority of cases, this tool will be used to create the new ROI centred on the "feature" from the previous Step.
 **/
class ViewCreatorCentredOnIParticleROITool: public extends<AthAlgTool, IViewCreatorROITool>
{
public:
  ViewCreatorCentredOnIParticleROITool(const std::string& type, const std::string& name, const IInterface* parent);

  virtual ~ViewCreatorCentredOnIParticleROITool() = default;

   virtual StatusCode initialize() override;

  /**
   * @brief Tool interface method.
   **/
  virtual StatusCode attachROILinks(TrigCompositeUtils::DecisionContainer& decisions, const EventContext& ctx) const override;

  SG::WriteHandleKey< TrigRoiDescriptorCollection > m_roisWriteHandleKey {this,"RoisWriteHandleKey","",
    "Name of the ROI collection produced by this tool."};

  Gaudi::Property< std::string > m_iParticleLinkName{this,"IParticleLinkName","feature",
    "Name of linked IParticle object to centre the new ROI on. Normally the 'feature' from the previous Step."};

  Gaudi::Property< double > m_roiEtaWidth{this,"RoIEtaWidth",0.1,
    "Extent of the ROI in eta from its centre"};

  Gaudi::Property< double > m_roiPhiWidth{this,"RoIPhiWidth",0.1,
    "Extent of the ROI in phi from its centre"};

  Gaudi::Property< double > m_roiZedWidth{this,"RoIZedWidth",-999,
    "Extent of the ROI in zed from its centre"};

  Gaudi::Property< double > m_roiZedSinThetaFlag{this,"RoIZedSinThetaFlag", false,
    "use a delta zed/sin\theta width rather than just a zed width"};

};

#endif //> !DECISIONHANDLING_VIEWCREATORCENTREDONIPARTICLEROITOOL_H
