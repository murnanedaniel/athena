/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#ifndef ALFA_RAWDATAPROVIDERTOOL_CHARGE_H
#define ALFA_RAWDATAPROVIDERTOOL_CHARGE_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"

#include "ByteStreamData/RawEvent.h" 

#include "ALFA_RawEv/ALFA_RawDataContainer_charge.h"
#include "ALFA_RawEv/ALFA_RawDataCollection_charge.h"
#include "ALFA_RawEv/ALFA_RawData_charge.h"

#include "ALFA_Decoder_charge.h"

#include <stdint.h>

#include <vector>
#include <set>
#include <string>


class ALFA_RawData_charge;
class ALFA_RawDataCollection_charge;
class ALFA_RawDataContainer_charge;


// the tool to decode a ROB fragment

class ALFA_RawDataProviderTool_charge : public AthAlgTool
{

 public:
   
  //! constructor
  ALFA_RawDataProviderTool_charge(const std::string& type, const std::string& name, const IInterface* parent);

  //! destructor
  virtual ~ALFA_RawDataProviderTool_charge();

  //! initialize
  virtual StatusCode initialize() override;

   //! this is the main decoding method
  StatusCode convert_charge(std::vector<const OFFLINE_FRAGMENTS_NAMESPACE::ROBFragment*>& vecRobs,ALFA_RawDataContainer_charge* rdoCont);
  
private:

  ToolHandle<ALFA_Decoder_charge>  m_decoder_charge{this, "Decoder_charge", "ALFA_Decoder_charge"};

};

#endif

