/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack


#ifndef SELECTION_HELPERS__OUT_OF_VALIDITY_HELPER_H
#define SELECTION_HELPERS__OUT_OF_VALIDITY_HELPER_H

#include <AsgMessaging/AsgMessagingForward.h>
#include <AthContainers/AuxElement.h>
#include <CxxUtils/AthUnlikelyMacros.h>
#include <SelectionHelpers/ISelectionWriteAccessor.h>
#include <xAODBase/IParticle.h>
#include <memory>

class StatusCode;

namespace CP
{
  class CorrectionCode;


  /// \brief the action to perform on encountering an
  /// OutOfValidityRange in \ref OutOfValidityHelper
  ///
  /// This is in addition to possibly setting a selection decoration
  /// (if that has been configured.
  ///
  /// This is not a member enum of \ref OutOfValidityHelper to void
  /// problems with member enums in dictionary generation, etc.

  enum class OutOfValidityAction
  {
    /// \brief print an error message and return a failure status code
    /// (triggering an abort)
    ///
    /// This is currently the default, as it forces users to think
    /// about what they want to do with OutOfValidityRange results.
    /// That is fairly safe, as most tools don't actually report
    /// OutOfValidityRange and we don't want to add a lot of
    /// meaningless selection decorations to each object.
    ABORT,

    /// \brief print a warning message and return a success status
    /// code.
    ///
    /// This should only be used if OutOfValidity is a *rare* event,
    /// or otherwise it will completely clobber the log file.
    WARNING,

    /// \brief don't print anything and return success
    ///
    /// This should (normally) be combined with a selection decoration
    /// that records OutOfValidity results, allowing to retrieve that
    /// information subsequently.
    SILENT
  };



  /// \brief a helper to translate a \ref CP::CorrectionCode into a
  /// \ref ::StatusCode
  ///
  /// The prolem is OutOfValidityRange which does not have an
  /// equivalent in StatusCode and which does not have a unique,
  /// correct handling in all situations.  This helper allows to
  /// configure a variety of behaviors via properties.

  class OutOfValidityHelper final : public asg::AsgMessagingForward
  {
    /// \brief standard constructor
  public:
    template<typename T>
    OutOfValidityHelper (T *owner, const std::string& propertyName = "outOfValidity",
                         const std::string& propertyDescription = "how to handle out of validity results");


    /// \brief standard initialize
  public:
    ::StatusCode initialize ();

    /// \brief check the correction code and do the proper thing
  public:
    ::StatusCode check (const xAOD::IParticle& particle,
                        const CP::CorrectionCode& code,
                        const char *context) const;


    /// \brief the action to take
  private:
    unsigned m_action {unsigned (OutOfValidityAction::ABORT)};

    /// \brief the accessor if we apply one
  private:
    std::unique_ptr<ISelectionWriteAccessor> m_accessor;

    /// \brief the decoration name we use (if we have one)
  private:
    std::string m_decorationName;

    /// \brief whether we have been initialized
    ///
    /// This is only used in debug mode to indicate a programming
    /// fault.  Otherwise it is too easy for users to forget to
    /// initialize this object.
  private:
    bool m_isInitialized = false;
  };



  template<typename T> OutOfValidityHelper ::
  OutOfValidityHelper (T *owner, const std::string& propertyName,
                       const std::string& propertyDescription)
    : asg::AsgMessagingForward (owner)
  {
    owner->declareProperty (propertyName, m_action,
                            propertyDescription);
    owner->declareProperty (propertyName + "Deco", m_decorationName,
                            "decoration to set alongside action described by " + propertyName);
  }
}

/// \brief a helper check macro to work with \ref OutOfValidityHelper
#define ANA_CHECK_CORRECTION(helper,object,expr)                \
  { if (ATH_UNLIKELY((helper).check ((object), (expr), #expr).isFailure())) \
      return StatusCode::FAILURE; }

#endif
