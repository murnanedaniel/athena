/*
 *   Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
 *   */

//////////////////////////////////////////////////////////
/// class ICPJetCorrectionTool
///
/// Interface class for smearing the jet mass scale and resolution of large-R jets
/// It allows the user to derive the systematic uncertainties associated with the JMS and the JMR
///
/// For information, see the Twiki:
/// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/FFJetSmearingTool
///
//////////////////////////////////////////////////////////

#ifndef JETCPINTERFACES_ICPJETCORRECTIONTOOL_H
#define JETCPINTERFACES_ICPJETCORRECTIONTOOL_H

#include "PATInterfaces/CorrectionCode.h"
#include "PATInterfaces/ISystematicsTool.h"

#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
#include "xAODEventInfo/EventInfo.h"

class ICPJetCorrectionTool : virtual public asg::IAsgTool,
                             virtual public CP::ISystematicsTool
{
    // Interface declaration
    ASG_TOOL_INTERFACE(ICPJetCorrectionTool)

    public:

        /// Virtual destructor
//        virtual ~ICPJetCorrectionTool()=default;

        /// Apply a systematic variation 
        virtual CP::CorrectionCode applyCorrection(xAOD::Jet& jet_reco) = 0;

}; // class ICPJetCorrectionTool

#endif /* JETINTERFACE_IJETCPCORRECTIONTOOL_H  */
