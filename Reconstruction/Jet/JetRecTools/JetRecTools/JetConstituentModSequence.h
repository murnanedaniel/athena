// this file is -*- C++ -*-

/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//JetConstituentModSequence.h

#ifndef JETRECTOOLS_JETCONSTITUENTMODSEQUENCE_H
#define JETRECTOOLS_JETCONSTITUENTMODSEQUENCE_H

//
// Michael Nelson, CERN & Univesity of Oxford
// February, 2016

#include <string>
#include "xAODBase/IParticleHelpers.h"
#include "xAODBase/IParticleContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODPFlow/TrackCaloClusterContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODPFlow/PFOContainer.h"
#include "xAODPFlow/FlowElementContainer.h"

#include "AsgTools/AsgTool.h"
#include "JetInterface/IJetExecuteTool.h"
#include "JetInterface/IJetConstituentModifier.h"
#include "AsgTools/ToolHandleArray.h"
#include "AsgDataHandles/ReadHandleKey.h"
#include "AsgDataHandles/WriteHandleKey.h"
#include "AsgDataHandles/ReadHandle.h"
#include "xAODCore/ShallowCopy.h"
#include "AsgTools/PropertyWrapper.h"

#ifndef XAOD_ANALYSIS
#include "AthenaMonitoringKernel/GenericMonitoringTool.h"
#endif

class JetConstituentModSequence: public asg::AsgTool, virtual public IJetExecuteTool {
  // Changed from IJetExecute
  ASG_TOOL_CLASS(JetConstituentModSequence, IJetExecuteTool)
  public:
  JetConstituentModSequence(const std::string &name); // MEN: constructor 
  StatusCode initialize();
  int execute() const;

protected:
  Gaudi::Property<std::string> m_inputContainer {this, "InputContainer", "", "The input container for the sequence"};
  Gaudi::Property<std::string> m_outputContainer = {this, "OutputContainer", "", "The output container for the sequence"};

  // P-A : the actual type
  // Define as a basic integer type because Gaudi
  // doesn't support arbitrary property types
  unsigned short m_inputType; // 
  
  
  ToolHandleArray<IJetConstituentModifier> m_modifiers{this , "Modifiers" , {} , "List of constit modifier tools."};

#ifndef XAOD_ANALYSIS
  ToolHandle<GenericMonitoringTool> m_monTool{this,"MonTool","","Monitoring tool"};
#endif
  
  bool m_saveAsShallow = true;

  // note: not all keys will be used for a particular instantiation
  SG::ReadHandleKey<xAOD::CaloClusterContainer> m_inClusterKey{this, "InClusterKey", "", "ReadHandleKey for unmodified CaloClusters"};
  SG::WriteHandleKey<xAOD::CaloClusterContainer> m_outClusterKey{this, "OutClusterKey", "", "WriteHandleKey for modified CaloClusters"};

  SG::ReadHandleKey<xAOD::TrackCaloClusterContainer> m_inTCCKey{this, "InTCCKey", "", "ReadHandleKey for unmodified TrackCaloClusters"};
  SG::WriteHandleKey<xAOD::TrackCaloClusterContainer> m_outTCCKey{this, "OutTCCKey", "", "WriteHandleKey for modified TrackCaloClusters"};

  SG::ReadHandleKey<xAOD::PFOContainer> m_inChargedPFOKey{this, "InChargedPFOKey", "", "ReadHandleKey for modified Charged PFlow Objects"};
  SG::WriteHandleKey<xAOD::PFOContainer> m_outChargedPFOKey{this, "OutChargedPFOKey", "", "WriteHandleKey for modified Charged PFlow Objects"};

  SG::ReadHandleKey<xAOD::PFOContainer> m_inNeutralPFOKey{this, "InNeutralPFOKey", "", "ReadHandleKey for modified Neutral PFlow Objects"};
  SG::WriteHandleKey<xAOD::PFOContainer> m_outNeutralPFOKey{this, "OutNeutralPFOKey", "", "WriteHandleKey for modified Neutral PFlow Objects"};

  SG::WriteHandleKey<xAOD::PFOContainer> m_outAllPFOKey{this, "OutAllPFOKey", "", "WriteHandleKey for all modified PFlow Objects"};

  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inChargedFEKey{this, "InChargedFEKey", "", "ReadHandleKey for modified Charged FlowElements"};
  SG::WriteHandleKey<xAOD::FlowElementContainer> m_outChargedFEKey{this, "OutChargedFEKey", "", "WriteHandleKey for modified Charged FlowElements"};

  SG::ReadHandleKey<xAOD::FlowElementContainer> m_inNeutralFEKey{this, "InNeutralFEKey", "", "ReadHandleKey for modified Neutral FlowElements"};
  SG::WriteHandleKey<xAOD::FlowElementContainer> m_outNeutralFEKey{this, "OutNeutralFEKey", "", "WriteHandleKey for modified Neutral FlowElements"};

  SG::WriteHandleKey<xAOD::FlowElementContainer> m_outAllFEKey{this, "OutAllFEKey", "", "WriteHandleKey for all modified FlowElements"};

  StatusCode copyModRecordPFO() const;
  StatusCode copyModRecordFE() const;

  /// helper function to cast, shallow copy and record a container.

  template<class T>
  StatusCode
  copyModRecord(const SG::ReadHandleKey<T>&,
                const SG::WriteHandleKey<T>&) const;
};

template<class T>
StatusCode
JetConstituentModSequence::copyModRecord(const SG::ReadHandleKey<T>& inKey,
                                         const SG::WriteHandleKey<T>& outKey) const {

  /* Read in a container of (type is template parameter),
     optionally modify the elements of this container, and store.
     This puts a (modified) copy of the container  into storegate.
  */
  
  auto inHandle = makeHandle(inKey);
  if(!inHandle.isValid()){
    ATH_MSG_WARNING("Unable to retrieve input container from " << inKey.key());
    return StatusCode::FAILURE;
  }

  std::pair< T*, xAOD::ShallowAuxContainer* > newconstit =
    xAOD::shallowCopyContainer(*inHandle);    
  newconstit.second->setShallowIO(m_saveAsShallow);

  for (auto t : m_modifiers) {ATH_CHECK(t->process(newconstit.first));}

  auto handle = makeHandle(outKey);
  ATH_CHECK(handle.record(std::unique_ptr<T>(newconstit.first),
                          std::unique_ptr<xAOD::ShallowAuxContainer>(newconstit.second)));
  
  xAOD::setOriginalObjectLink(*inHandle, *handle);
  
  return StatusCode::SUCCESS;
}

#endif
