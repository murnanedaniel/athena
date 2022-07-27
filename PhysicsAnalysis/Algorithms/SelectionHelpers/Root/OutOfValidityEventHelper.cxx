/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

/// @author Nils Krumnack


//
// includes
//

#include <SelectionHelpers/OutOfValidityEventHelper.h>

#include <AsgMessaging/MessageCheck.h>
#include <AsgMessaging/StatusCode.h>
#include <PATInterfaces/CorrectionCode.h>

//
// method implementations
//

namespace CP
{
  ::StatusCode OutOfValidityEventHelper ::
  initialize ()
  {
    m_isInitialized = true;
    return StatusCode::SUCCESS;
  }


  ::StatusCode OutOfValidityEventHelper ::
  check (const CP::CorrectionCode& code,
         const char *context) const
  {
    assert (m_isInitialized);

    switch (code)
    {
    case CP::CorrectionCode::Ok:
      return StatusCode::SUCCESS;
    case CP::CorrectionCode::Error:
      return StatusCode::FAILURE;
    case CP::CorrectionCode::OutOfValidityRange:
      switch (OutOfValidityAction (m_action))
      {
      case OutOfValidityAction::ABORT:
        ANA_MSG_ERROR ("encountered OutOfValidity: " << context);
        return StatusCode::FAILURE;
      case OutOfValidityAction::WARNING:
        ANA_MSG_WARNING ("encountered OutOfValidity: " << context);
        return StatusCode::SUCCESS;
      case OutOfValidityAction::SILENT:
        return StatusCode::SUCCESS;
      }
    }
    ANA_MSG_ERROR (__FILE__ << ":" << __LINE__ << ": invalid enum value encountered " << code << " " << int (m_action));
    return StatusCode::FAILURE;
  }

}
