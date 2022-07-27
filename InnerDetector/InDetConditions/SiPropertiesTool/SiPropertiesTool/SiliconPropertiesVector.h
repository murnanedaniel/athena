/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SIPROPERTIESTOOL_SILICONPROPERTIESVECTOR_H
#define SIPROPERTIESTOOL_SILICONPROPERTIESVECTOR_H

///////////////////////////////////////////////////////////////////
//   Header file for class SiliconPropertiesVector for SiPropertiesTool
///////////////////////////////////////////////////////////////////

#include <vector>
#include "SiPropertiesTool/SiliconProperties.h"
#include "Identifier/IdentifierHash.h"
#include "AthenaKernel/CLASS_DEF.h"

namespace InDet {

class SiliconPropertiesVector
{
public:
 
  SiliconPropertiesVector();
  ~SiliconPropertiesVector();

  void resize(IdentifierHash::value_type size);
  const InDet::SiliconProperties& getSiProperties(const IdentifierHash& elementHash) const;
  void setConditions(const IdentifierHash& elementHash, const double temperature, const double electricField);

private:
  std::vector<InDet::SiliconProperties> m_data;
  InDet::SiliconProperties m_defaultProperties;
};

} // namespace InDet

CLASS_DEF(InDet::SiliconPropertiesVector, 10146683 , 1)

#include "AthenaKernel/CondCont.h"
CONDCONT_MIXED_DEF(InDet::SiliconPropertiesVector, 246562685);

#endif // SIPROPERTIESTOOL_SILICONPROPERTIESVECTOR_H
