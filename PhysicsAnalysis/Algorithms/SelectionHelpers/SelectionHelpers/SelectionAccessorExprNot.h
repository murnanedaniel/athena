/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef SELECTION_HELPERS__SELECTION_ACCESSOR_EXPR_AND_H
#define SELECTION_HELPERS__SELECTION_ACCESSOR_EXPR_AND_H

#include "SelectionHelpers/SelectionAccessorExprBase.h"

#include <memory>

namespace CP {

/// @brief SelectionAccessor which implements a logical NOT
class SelectionAccessorExprNot : public SelectionAccessorExprBase {
 public:
  /// @param child The selection accessor that is to be negated
  SelectionAccessorExprNot(std::unique_ptr<ISelectionReadAccessor> child);

  /// @param element AuxElement to evaluate the selection on
  /// @return Result of the negated result of @c child
  virtual bool getBool(const SG::AuxElement &element,
                       const CP::SystematicSet *sys) const override;

  /// @return Returns a readable and parseable representation
  virtual std::string label() const override;

  virtual CP::SystematicSet
  getInputAffecting (const ISystematicsSvc& svc,
                     const std::string& objectName) const override;

  virtual StatusCode
  fillSystematics (const ISystematicsSvc& svc,
                   const std::vector<CP::SystematicSet>& sysList,
                   const std::string& objectName) override;

 private:
  std::unique_ptr<ISelectionReadAccessor> m_child;
};


}  // namespace CP

#endif
