/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack


#ifndef SELECTION_HELPERS__OUT_OF_VALIDITY_EVENT_HELPER_H
#define SELECTION_HELPERS__OUT_OF_VALIDITY_EVENT_HELPER_H

#include <AsgMessaging/AsgMessagingForward.h>
#include <AthContainers/AuxElement.h>
#include <CxxUtils/AthUnlikelyMacros.h>
#include <SelectionHelpers/ISelectionWriteAccessor.h>
#include <SelectionHelpers/OutOfValidityHelper.h>
#include <xAODBase/IParticle.h>
#include <memory>

class StatusCode;

namespace CP
{
  class CorrectionCode;



  /// \brief a helper to translate a \ref CP::CorrectionCode into a
  /// \ref ::StatusCode
  ///
  /// The prolem is OutOfValidityRange which does not have an
  /// equivalent in StatusCode and which does not have a unique,
  /// correct handling in all situations.  This helper allows to
  /// configure a variety of behaviors via properties.

  class OutOfValidityEventHelper final : public asg::AsgMessagingForward
  {
    /// \brief standard constructor
  public:
    template<typename T>
    OutOfValidityEventHelper (T *owner, const std::string& propertyName = "outOfValidity",
                         const std::string& propertyDescription = "how to handle out of validity results");


    /// \brief standard initialize
  public:
    ::StatusCode initialize ();

    /// \brief check the correction code and do the proper thing
  public:
    ::StatusCode check (const CP::CorrectionCode& code,
                        const char *context) const;


    /// \brief the action to take
  private:
    unsigned m_action {unsigned (OutOfValidityAction::ABORT)};

    /// \brief whether we have been initialized
    ///
    /// This is only used in debug mode to indicate a programming
    /// fault.  Otherwise it is too easy for users to forget to
    /// initialize this object.
  private:
    bool m_isInitialized = false;
  };



  template<typename T> OutOfValidityEventHelper ::
  OutOfValidityEventHelper (T *owner, const std::string& propertyName,
                       const std::string& propertyDescription)
    : asg::AsgMessagingForward (owner)
  {
    owner->declareProperty (propertyName, m_action,
                            propertyDescription);
  }
}

/// \brief a helper check macro to work with \ref OutOfValidityEventHelper
#define ANA_CHECK_CORRECTION_EVENT(helper,expr)       \
  { if (ATH_UNLIKELY((helper).check ((expr), #expr).isFailure())) \
      return StatusCode::FAILURE; }

#endif
