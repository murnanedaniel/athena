// -*- C++ -*-

/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file ISCT_ByteStreamErrorsTool.h
 * interface file for tool that keeps track of errors in the bytestream.
 * @author Susumu.Oda@cern.ch
 **/

#ifndef ISCT_ByteStreamErrorsTool_h
#define ISCT_ByteStreamErrorsTool_h

#include "InDetConditionsSummaryService/InDetHierarchy.h"
#include "SCT_ConditionsTools/ISCT_ConditionsTool.h"
#include "SCT_ConditionsData/SCT_ByteStreamErrors.h"
#include "InDetByteStreamErrors/IDCInDetBSErrContainer.h"
#include "GaudiKernel/EventContext.h"

#include <set>

class Identifier;
class IdentifierHash;

/**
 * @class SCT_ByteStreamErrorsTool
 * Tool that keeps track of modules that give rise to errors in the bytestram.
 **/
namespace SCT_ByteStreamErrors {
  //!< @brief helper to be used in clients to fetch error information
  inline bool hasError(IDCInDetBSErrContainer::ErrorCode errWord,  errorTypes errType ) { return errWord & (1<<errType); }

  //!< @brief helper to set the error: @example errors[hashId] = addError( errors[hashId], PixelByteStreamErrors::Invalid )
  inline void addError(IDCInDetBSErrContainer::ErrorCode& errWord,  errorTypes errType ) { errWord |= (1<<errType); }

  //!< @brief for cases when error doe snot need to be accumulated
  inline IDCInDetBSErrContainer::ErrorCode makeError( errorTypes errType ) { return (1<<errType); }
}


class ISCT_ByteStreamErrorsTool: virtual public ISCT_ConditionsTool {

 public:
  //@name Tool methods
  //@{

  virtual ~ISCT_ByteStreamErrorsTool() = default;

  /// Creates the InterfaceID and interfaceID() method
  DeclareInterfaceID(ISCT_ByteStreamErrorsTool, 1, 0);
  //@}
  
  virtual std::set<IdentifierHash> getErrorSet(int errorType) const =0;
  virtual std::set<IdentifierHash> getErrorSet(int errorType, const EventContext& ctx) const =0;

  /** Temporary status of chips for a particular module (packed as 1st 12 bits of unsigned int) */
  virtual unsigned int tempMaskedChips(const Identifier& moduleId) const =0;
  virtual unsigned int tempMaskedChips(const Identifier& moduleId, const EventContext& ctx) const =0;
  /** Status ABCD errors of chips for a particular module (packed as 1st 12 bits of unsigned int) */
  virtual unsigned int abcdErrorChips(const Identifier& moduleId) const =0;
  virtual unsigned int abcdErrorChips(const Identifier& moduleId, const EventContext& ctx) const =0;

 private:

};

#endif // ISCT_ByteStreamErrorsTool_h
