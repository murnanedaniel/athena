/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/


#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Service.h"

#include "CaloIdentifier/CaloCell_ID.h"

#include "CaloUtils/CaloTowerBuilderTool.h"

#include "LArRecUtils/LArTowerBuilderTool.h"

#include <string>

LArTowerBuilderTool::LArTowerBuilderTool(const std::string& name,
					 const std::string& type,
					 const IInterface* parent)
  : CaloTowerBuilderTool(name,type,parent)
{ }

LArTowerBuilderTool::~LArTowerBuilderTool()
{ }

StatusCode LArTowerBuilderTool::initializeTool()
{
  // services
  MsgStream log(msgSvc(),name());

  // allow only LAREM and LARHEC cells!
  for ( size_t iCalos=0; iCalos<m_includedCalos.size(); iCalos++ )
    {
      if ( m_includedCalos[iCalos] == "LAREM" )
	{
	  m_caloIndices.push_back(CaloCell_ID::LAREM);
	}
      else if ( m_includedCalos[iCalos] == "LARHEC" )
	{
	  m_caloIndices.push_back(CaloCell_ID::LARHEC);
	}
      else if ( m_includedCalos[iCalos] == "LARFCAL" )
	{
	  log << MSG::INFO
	      << "use LArFCalTowerBuilderTool for the FCal - request ignored"
	      << endreq;
	}
    }

  // check setup
  return this->checkSetup(log);
}
