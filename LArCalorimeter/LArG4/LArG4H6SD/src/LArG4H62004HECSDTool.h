/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef LARG4H62004HECSDTOOL_H
#define LARG4H62004HECSDTOOL_H

#include "LArG4Code/LArG4SDTool.h"
#include <string>
#include <vector>

#include "StoreGate/WriteHandle.h"
#include "LArSimEvent/LArHitContainer.h"

class LArG4SimpleSD;

class LArG4H62004HECSDTool : public LArG4SDTool
{
 public:
  // Constructor
  LArG4H62004HECSDTool(const std::string& type, const std::string& name, const IInterface *parent);
    
  // Destructor
  virtual ~LArG4H62004HECSDTool() {}

  // Method in which all the SDs are created and assigned to the relevant volumes
  StatusCode initializeSD() override final;

  // Calls down to all the SDs to get them to pack their hits into a central collection
  StatusCode Gather() override final;

 private:
  // The actual hit container - here because the base class is for both calib and standard SD tools
  SG::WriteHandle<LArHitContainer> m_HitColl;

  LArG4SimpleSD* m_wheelSD;

};

#endif
