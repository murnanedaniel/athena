/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "G4TruthStrategies/LLPStrategy.h"

#include "MCTruthBase/TruthStrategyUtils.h"
//#include "SimHelpers/StepHelper.h"
//#include "MCTruth/EventInformation.h"

#include "G4Step.hh"
#include "G4ProcessType.hh"
#include "G4VProcess.hh"

#include <iostream>

static LLPStrategy Strategy("LLP");

LLPStrategy::LLPStrategy(const std::string& n)
  : TruthStrategy(n)
{;}

bool LLPStrategy::AnalyzeVertex(const G4Step* aStep)
{
  // Thread-safety issue
  //static StepHelper sHelper;
  //sHelper.SetStep(aStep);
  //G4int pSubType = sHelper.GetProcessSubType();

  G4int pSubType =
    aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessSubType();

  // All kinds of decay processes are included here...
  if(pSubType > 200 && pSubType < 299) {

    // Check if this is a sparticle - if not, return
    int pdg = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    if ( !isSUSYParticle( abs(pdg) ) ) return false;

    // Get the list of secondaries in the current step
    const auto& tracks = *aStep->GetSecondaryInCurrentStep();

    // Save the vertex using TruthStrategyUtils
    TruthStrategyUtils::saveSecondaryVertex(nullptr, aStep->GetPostStepPoint(), tracks);

    // TODO: Remove these uncontrolled printouts!
    std::cout << "ACHLLP: saved a truth vertex for pdg " << pdg << std::endl;
    return true;
  }

  // This next bit could easily be a separate truth strategy

  G4ProcessType pType =
    aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessType();

  // Save all interactions for user-defined processes, like rhadron interactions
  if(fUserDefined == pType) {

    // Get the list of secondaries in the current step
    const auto& tracks = *aStep->GetSecondaryInCurrentStep();

    // Save the vertex using TruthStrategyUtils
    G4TrackStatus tStatus = aStep->GetTrack()->GetTrackStatus();
    G4Track* primTrack = (tStatus==fAlive) ? aStep->GetTrack() : nullptr;
    TruthStrategyUtils::saveSecondaryVertex(primTrack, aStep->GetPostStepPoint(), tracks);

    // TODO: Remove these uncontrolled printouts!
    int pdg = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    std::cout << "ACHLLP: saved a truth interaction fUserDefined for pdg " << pdg << std::endl;

    return true;
  }

  return false;	
}

bool LLPStrategy::isSUSYParticle(const int id) const
{
  return ( (id>1000000 && id<1000007) || (id>1000010 && id<1000017) ||
           (id>2000000 && id<2000007) || (id>2000010 && id<2000017) || // squarks and sleptons
           (id>1000020 && id<1000026) || (id>1000034 && id<1000040) || // higgses etc
           id==1000512 || id==1000522 || id==1000991 || id==1000993 || // R-hadrons
           id==1000612 || id==1000622 || id==1000632 || id==1000642 || id==1000652 || id==1005211 ||
           id==1006113 || id==1006211 || id==1006213 || id==1006223 || id==1006311 ||
           id==1006313 || id==1006321 || id==1006323 || id==1006333 ||
           id==1009111 || id==1009113 || id==1009211 || id==1009213 || id==1009311 ||
           id==1009313 || id==1009321 || id==1009323 || id==1009223 || id==1009333 ||
           id==1092112 || id==1091114 || id==1092114 || id==1092212 || id==1092214 || id==1092224 ||
           id==1093114 || id==1093122 || id==1093214 || id==1093224 || id==1093314 || id==1093324 || id==1093334 );
}

