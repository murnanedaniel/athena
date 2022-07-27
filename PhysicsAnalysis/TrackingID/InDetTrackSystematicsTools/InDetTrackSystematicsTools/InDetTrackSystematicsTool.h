// -*- c++ -*-
/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef INDETTRACKSYSTEMATICSTOOLS_INDETTRACKSYSTEMATICSTOOL_H
#define INDETTRACKSYSTEMATICSTOOLS_INDETTRACKSYSTEMATICSTOOL_H

#include <string>
#include <map>
#include <memory>

#include "InDetTrackSystematicsTools/InDetTrackSystematics.h"
#include "AsgTools/AsgTool.h"
#include "PATInterfaces/ISystematicsTool.h"
#include "PATInterfaces/SystematicSet.h"

#include <TFile.h>

namespace InDet {

  class InDetTrackSystematicsTool : public virtual CP::ISystematicsTool, public asg::AsgTool {

    ASG_TOOL_CLASS2( InDetTrackSystematicsTool,
                     CP::ISystematicsTool, CP::IReentrantSystematicsTool )
  public:
    InDetTrackSystematicsTool( const std::string& );
    virtual ~InDetTrackSystematicsTool() = default;

    virtual StatusCode initialize() override;

    /// returns: whether the tool is affected by the systematic
    virtual bool isAffectedBySystematic( const CP::SystematicVariation& ) const override;
    /// returns: list of systematics this tool can be affected by
    virtual CP::SystematicSet affectingSystematics() const override = 0;
    /// returns: list of recommended systematics to use with this tool
    virtual CP::SystematicSet recommendedSystematics() const override;
    /// configure the tool to apply a given list of systematic variations
    virtual StatusCode applySystematicVariation( const CP::SystematicSet& ) override;


  protected:
    /// open and return a file with the given name.
    std::unique_ptr<TFile> getFile( const std::string& ) const;

    /// a function to initialize an object from a root file
    template <class T> StatusCode initObject(T*& obj, std::string rootFileName, std::string objName) const;

    // a map from a general set to a set that is filtered for the ones we use
    // note that SystematicSet caches its hashes so this is probably not as slow as it might seem
    std::unordered_map< CP::SystematicSet, CP::SystematicSet > m_sysFilterMap;

    const CP::SystematicSet* m_activeSysts = nullptr;

    bool isActive( TrackSystematic ) const;

  };

}

// must define template function in header
template <class T>
StatusCode InDet::InDetTrackSystematicsTool::initObject(T*& obj, std::string rootFileName, std::string objName) const
{
  if (obj != nullptr) ATH_MSG_WARNING( obj->GetName() << " is not null, yet we are now attempting to initialize from " << rootFileName );
  std::unique_ptr<TFile> F = getFile(rootFileName);
  if(!F || F->IsZombie()) {
    ATH_MSG_ERROR( "Could not open file " << rootFileName );
    return StatusCode::FAILURE;
  }
  T* tempObj = nullptr;
  F->GetObject(objName.data(), tempObj);
  if(tempObj==nullptr) {
    ATH_MSG_ERROR( "Could not retrieve " << objName << " from file " << rootFileName );
    return StatusCode::FAILURE;
  }
  obj = static_cast<T*>(tempObj->Clone());
  obj->SetDirectory(0);
  F->Clear();
  F->Close();
  return StatusCode::SUCCESS;
}

#endif
