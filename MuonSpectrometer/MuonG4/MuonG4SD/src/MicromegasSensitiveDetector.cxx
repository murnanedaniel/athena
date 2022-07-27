/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

#include "MicromegasSensitiveDetector.h"
#include "MuonSimEvent/MicromegasHitIdHelper.h"
#include "MCTruth/TrackHelper.h"
#include "G4Geantino.hh"
#include "G4ChargedGeantino.hh"

#include "G4Track.hh"

#include "GeoPrimitives/CLHEPtoEigenConverter.h"

#include <string>

// construction/destruction
MicromegasSensitiveDetector::MicromegasSensitiveDetector(const std::string& name, const std::string& hitCollectionName)
  : G4VSensitiveDetector( name )
  , m_MMSimHitCollection( hitCollectionName )
{
  m_muonHelper = MicromegasHitIdHelper::GetHelper();
  //m_muonHelper->PrintFields();
}

// Implemenation of memebr functions
void MicromegasSensitiveDetector::Initialize(G4HCofThisEvent*) 
{
  if (!m_MMSimHitCollection.isValid()) m_MMSimHitCollection = std::make_unique<MMSimHitCollection>();
}

G4bool MicromegasSensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory* /*ROHist*/) 
{
  G4Track* currentTrack = aStep->GetTrack();
  int charge=currentTrack->GetDefinition()->GetPDGCharge();
  bool geantinoHit = (currentTrack->GetDefinition()==G4Geantino::GeantinoDefinition()) ||
                     (currentTrack->GetDefinition()==G4ChargedGeantino::ChargedGeantinoDefinition());

  if (!charge && (!geantinoHit)) return false;
  //  G4cout << "\t\t MicromegasSD: Hit in a sensitive layer!!!!! " << G4endl;
  G4StepPoint* postStep=aStep->GetPostStepPoint();
  const G4Step* post_Step=aStep->GetTrack()->GetStep();

  Amg::Vector3D position = Amg::Hep3VectorToEigen( postStep->GetPosition() );
  
  int pdgCode=currentTrack->GetDefinition()->GetPDGEncoding();

  float globalTime=postStep->GetGlobalTime();
  float eKin=postStep->GetKineticEnergy();
  if (eKin<= 0. && (!geantinoHit)) return false;

  Amg::Vector3D direction = Amg::Hep3VectorToEigen( postStep->GetMomentumDirection() );
  float depositEnergy=post_Step->GetTotalEnergyDeposit();
  
  if (depositEnergy<0.0001 && (!geantinoHit)) return false;
  
  const G4TouchableHistory* touchHist = static_cast<const G4TouchableHistory*>(aStep->GetPreStepPoint()->GetTouchable());
  
  // int iDepth=touchHist->GetHistoryDepth();
  //  G4cout << "\t\t\t\t Touchable history dump " << G4endl;
  int nLayer=touchHist->GetVolume(0)->GetCopyNo();
  std::string chName=touchHist->GetVolume(1)->GetLogicalVolume()->GetName();
  std::string subType=chName.substr(chName.find('-')+1);
  if (subType[0]!='M') G4cout << " something is wrong, this is no Micromegas!" << G4endl;
  std::string temp(&subType[1]);
  std::istringstream is(temp);
  int iRing;
  is>>iRing;
  // identifiers have eta naming 0-1, eta encoded in subtype is 1-2
  iRing--;
  // double phiDiff=2*M_PI;

  G4ThreeVector posH=postStep->GetPosition(); //posH is equivalent to position - eigen not used to avoid additional dependence on EventPrimitives
  if (subType[2]=='L') posH.rotateZ(M_PI/8.);
  double phiHit=posH.phi();
  if(phiHit<=0) phiHit+=2.*M_PI;
  int iPhi=1+int(phiHit/(M_PI/4.));
  iPhi*=2;
  if (subType[2]=='L') iPhi-=1;
  
  int iSide=1;
  if (position.z()<0) iSide=-1;
  
  int mLayer= atoi((subType.substr(3,1)).c_str());
  if (mLayer != 1 && mLayer !=2) G4cout << " something is wrong - multilayer index is " << mLayer << G4endl;
  
  //  G4cout << "\t\t Chamber "<<chName<<" subType "<<subType<<" layer nr. "<<nLayer<<" ring "<<iRing<<" sector "<<iPhi<<" side "<<iSide << G4endl;
  int MmId = m_muonHelper->BuildMicromegasHitId(subType, iPhi, iRing, mLayer,nLayer, iSide);
 
  TrackHelper trHelp(aStep->GetTrack());

  m_MMSimHitCollection->Emplace(MmId, globalTime,position,pdgCode,eKin,direction,depositEnergy,trHelp.GetParticleLink());

  //    G4cout << "MMs "<<m_muonHelper->GetStationName(MmId)
  // 	            << " "<<m_muonHelper->GetFieldValue("PhiSector")
  // 	            << " "<<m_muonHelper->GetFieldValue("ZSector")
  // 	            << " "<<m_muonHelper->GetFieldValue("MultiLayer")
  // 	            << " "<<m_muonHelper->GetFieldValue("Layer")
  // 	            << " "<<m_muonHelper->GetFieldValue("Side") << G4endl;  	    

  //G4cout << m_muonHelper->GetStationName(MmId)<<" "<<aHit->print() << G4endl;
  //  G4cout << aHit->print() << G4endl;

  return true;
}

