/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef MUONCONDTOOL_MDT_MAPCONVERSION_H
#define MUONCONDTOOL_MDT_MAPCONVERSION_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "MuonCondInterface/IMDT_MapConversion.h"
#include "MuonIdHelpers/MdtIdHelper.h"

class IIOVSvc;

class MDT_MapConversion : public AthAlgTool, virtual public IMDT_MapConversion {
public:
    MDT_MapConversion(const std::string& type, const std::string& name, const IInterface* parent);

    virtual StatusCode initialize();

    virtual const Identifier& ConvertToOffline(std::string_view OnlineId, bool quiet = false) const;

    //  const std::string OnlineName(Identifier OfflineId);
    // const Identifier OfflineName(std::string OnlineId);

private:
    std::map<std::string, Identifier> m_Chamber_Map;

    const MdtIdHelper* m_mdtIdHelper;

    std::string m_OnlineName;
    Identifier m_OfflineName;
    Identifier m_Online_empty;

    std::string m_chrono1;
    std::string m_chrono2;
    std::string m_chrono3;
    std::string m_chrono4;
};

#endif
