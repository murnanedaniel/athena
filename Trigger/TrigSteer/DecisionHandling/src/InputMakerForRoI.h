/*
  Copyright (C) 2002-2018 CERN for the benefit of the ATLAS collaboration
*/
#ifndef TRIGUPGRADETEST_INPUTMAKERFORROI_H
#define TRIGUPGRADETEST_INPUTMAKERFORROI_H 


#include <string>
#include "DecisionHandling/InputMakerBase.h"
#include "DecisionHandling/TrigCompositeUtils.h"
#include "AthContainers/ConstDataVector.h"
#include "StoreGate/ReadHandleKeyArray.h"
#include "TrigSteeringEvent/TrigRoiDescriptorCollection.h"

  /**
   * @class  InputMakerForRoI
   * @brief Used at the start of a sequence: retrieves filtered collection via menu decision from previous step and writes the RoI collection out directly so it can be used as input by the reco alg that follows in sequence.
   **/

  using namespace TrigCompositeUtils;
  
  class  InputMakerForRoI    : public ::InputMakerBase  { 
  public: 
    InputMakerForRoI( const std::string& name, ISvcLocator* pSvcLocator );
    virtual ~ InputMakerForRoI(); 
    virtual StatusCode  initialize() override;
    virtual StatusCode  execute_r(const EventContext&) const override;
    virtual StatusCode  finalize() override;

  private: 
     InputMakerForRoI();
 
    // Use this to customise the example for a different input feature
    typedef TrigRoiDescriptor FeatureOBJ;
    typedef TrigRoiDescriptorCollection FeatureContainer;

    SG::WriteHandleKey<TrigRoiDescriptorCollection> m_RoIs {this,"RoIs", "Unspecified", "Nam eof the RoIs extracted from the decisions"};

    StringProperty m_linkName   {this, "LinkName", "initialRoI",  "name of the link to the features in the decision, e.g. 'feature', 'initialRoI'"};


  }; 


#endif //> !TRIGUPGRADETEST_INPUTMAKERFORROI_H
