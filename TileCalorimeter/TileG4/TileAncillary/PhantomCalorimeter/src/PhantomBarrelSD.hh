/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

//************************************************************
//
// Class PhantomBarrelSD
// Sensitive detector for the phantom calorimeter in combined 2004
//
//************************************************************

#ifndef PhantomCalorimeter_PhantomBarrelSD_H
#define PhantomCalorimeter_PhantomBarrelSD_H

// Base class
#include "G4VSensitiveDetector.hh"

// use of the hits
#include "TileSimEvent/TileHitVector.h"
#include "StoreGate/WriteHandle.h"

// STL header
#include <string>

class G4Step;
class G4TouchableHistory;

class G4Step;
class G4HCofThisEvent;
class TileTBID;

class PhantomBarrelSD : public G4VSensitiveDetector
{
public:
  PhantomBarrelSD(const std::string& name, const std::string& hitCollectionName);
  ~PhantomBarrelSD();

  // Called from PhantomBarrelSDTool::SetupEvent
  void StartOfAthenaEvent ();
  void Initialize(G4HCofThisEvent*) override final;
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) override final;
  void EndOfAthenaEvent();

private:
  const TileTBID* m_tileTBID;

  static const int nCell = 8;

  int m_nhits[nCell];
  TileSimHit* m_hit[nCell];
  Identifier m_id[nCell];
  // The hits collections
  SG::WriteHandle<TileHitVector> m_HitColl;
};

#endif
