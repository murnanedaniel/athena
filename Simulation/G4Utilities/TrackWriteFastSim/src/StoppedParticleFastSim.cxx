/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "TrackWriteFastSim/StoppedParticleFastSim.h"
#include "G4FastSimulation/FastSimModelProxy.h"
#include "FadsSensitiveDetector/SensitiveDetectorCatalog.h"
#include "FadsSensitiveDetector/FadsSensitiveDetector.h"

#include "TrackWriteFastSim/TrackFastSimSD.h"

#include "G4FastTrack.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"

#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <iostream>

static FastSimModelProxy<StoppedParticleFastSim> stoppedparticlefastsim("StoppedParticleFastSim");

G4bool StoppedParticleFastSim::IsApplicable(const G4ParticleDefinition&)
{
  return true;
}

G4bool StoppedParticleFastSim::ModelTrigger(const G4FastTrack& fastTrack)
{
  // Trigger if the energy is below our threshold or if the time is over 150 ns
  int id = fastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();
  if (id<1000000 || id>1100000) return true;
  else if (isSUSYParticle(id)){
    G4Material * mat = fastTrack.GetPrimaryTrack()->GetMaterial();
    double minA=1500000.;
    for (unsigned int i=0;i<mat->GetNumberOfElements();++i){
      if (mat->GetElement(i) &&
          minA>mat->GetElement(i)->GetN()){
        minA=mat->GetElement(i)->GetN();
        //std::cout << "Setting min A !" << std::endl;
      } //else if (! mat->GetElement(i) ) std::cout << "Had a problem with " << mat->GetElement(i) << std::endl;
      //else std::cout << "Got an A from the material of: " << mat->GetElement(i)->GetN() << std::endl;
    }
    //std::cout << "Material with " << mat->GetNumberOfElements() << " elements... and vector size " << mat->GetElementVector()->size() << std::endl;
    //std::cout << "Looking at SUSY particle " << id << " with beta " << fastTrack.GetPrimaryTrack()->GetVelocity() << " and minA " << minA << " for " << 0.15*std::pow(minA,-2./3.) << std::endl;
    if (fastTrack.GetPrimaryTrack()->GetVelocity()<0.15*std::pow(minA,-2./3.)*CLHEP::c_light) return true;
    return false;
  }
  return true;
}

void StoppedParticleFastSim::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep)
{
  if (!m_init){
    m_init = true;

    FADS::SensitiveDetectorCatalog * fsdc = FADS::SensitiveDetectorCatalog::GetSensitiveDetectorCatalog();
    if (!fsdc) { 
      std::cout << "ERROR: StoppedParticleFastSim could not get sensitive detector catalog." << std::endl;
    } else {
      FADS::FadsSensitiveDetector * fsd = fsdc->GetSensitiveDetector("TrackFastSimSD");
      if (!fsd) { 
        std::cout << "StoppedParticleFastSim could not get TrackFastSimSD sensitive detector." << std::endl; 
      } else {
        m_fsSD = dynamic_cast<TrackFastSimSD*>(fsd);
        if (!m_fsSD) {
          std::cout << "StoppedParticleFastSim could not cast the SD." << std::endl; 
        } else { // succeeded in cast
          m_fsSD->SetCollectionName("StoppingPositions");
        }
      } // found the SD
    } // got the catalog
  }
  
  if (isSUSYParticle(fastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetDefinition()->GetPDGEncoding()) &&
      m_fsSD) {
    m_fsSD->WriteTrack( fastTrack.GetPrimaryTrack() , false , true );
  }
  fastStep.KillPrimaryTrack();
}

bool StoppedParticleFastSim::isSUSYParticle(const int id) const
{
  if (id==1000021 || id==1000005 || id==1000006 || id==1000512 || id==1000522 || id==1000991 || id==1000993 ||
      id==1000612 || id==1000622 || id==1000632 || id==1000642 || id==1000652 || id==1005211 ||
      id==1006113 || id==1006211 || id==1006213 || id==1006223 || id==1006311 ||
      id==1006313 || id==1006321 || id==1006323 || id==1006333 ||
      id==1009111 || id==1009113 || id==1009211 || id==1009213 || id==1009311 ||
      id==1009313 || id==1009321 || id==1009323 || id==1009223 || id==1009333 ||
      id==1092112 || id==1091114 || id==1092114 || id==1092212 || id==1092214 || id==1092224 ||
      id==1093114 || id==1093122 || id==1093214 || id==1093224 || id==1093314 || id==1093324 || id==1093334) 
    return true;
  return false;
}

